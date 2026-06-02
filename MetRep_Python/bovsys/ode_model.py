"""Core 22-state metabolic-reproductive ODE system.

This module intentionally excludes Dexa PK/PD. It translates the non-Dexa
behavior of ``BovSys_Equa_dexa_v3`` from MATLAB for states 1..22.
"""

from __future__ import annotations

import numpy as np


def _interp(t: float, at: np.ndarray, values: np.ndarray) -> float:
    """Linear interpolation with MATLAB-like endpoint extrapolation by clamp."""

    return float(np.interp(t, at, values))


def bovsys_rhs(
    t: float,
    y: np.ndarray,
    params: dict[str, float],
    at: np.ndarray,
    dmi: np.ndarray,
    milk: np.ndarray,
    mode: str = "non_lactating",
) -> np.ndarray:
    """Return dy/dt for the core BovSys model.

    Parameters
    ----------
    t:
        Simulation time in days.
    y:
        State vector with 22 states in ``initial_conditions.STATE_NAMES`` order.
    params:
        Parameter dictionary from ``default_parameters``.
    at, dmi, milk:
        Time grid and forcing vectors for dry matter intake and milk yield.
    mode:
        ``"lactating"`` enables the lactation OXT input term from MATLAB.
    """

    y = np.maximum(np.asarray(y, dtype=float), 0.0)
    f = np.zeros(22, dtype=float)

    (
        y_gnrhh,
        y_gnrh,
        y_fshp,
        y_fsh,
        y_lhp,
        y_lh,
        y_fol,
        y_pg,
        y_cl,
        y_p4,
        y_e2,
        y_ih,
        y_enz,
        y_ot1,
        y_iof,
        y_igf_b,
        y_ins,
        y_glu,
        y_fat,
        y_lv,
        y_glu_st,
        y_gluca,
    ) = y

    DMI = _interp(t, at, dmi)
    Milk = _interp(t, at, milk)

    # GnRH
    hm_e2_1 = params["hm_e2_1_threshold"] ** 2 / (
        params["hm_e2_1_threshold"] ** 2 + y_e2**2
    )
    hm_p4_1 = params["hm_p4_1_threshold"] ** 2 / (
        params["hm_p4_1_threshold"] ** 2 + y_p4**2
    )
    hm_p4_2 = params["hm_p4_2_scale"] * params["hm_p4_2_threshold"] ** 2 / (
        params["hm_p4_2_threshold"] ** 2 + y_p4**2
    )
    hp_e2_1 = params["hp_e2_1_scale"] * y_e2**5 / (
        params["hp_e2_1_threshold"] ** 5 + y_e2**5
    )
    gnrh_syn = params["gnrh_syn"] * (1.0 - y_gnrhh / params["gnrh_hypo_max"])
    gnrh_rel = (
        params["gnrh_release_scale"] * (hm_p4_1 + hm_e2_1 - hm_e2_1 * hm_p4_1)
        + hm_p4_2
    )
    f[0] = gnrh_syn - gnrh_rel * y_gnrhh
    f[1] = gnrh_rel * y_gnrhh * hp_e2_1 - params["gnrh_clearance"] * y_gnrh

    # FSH with insulin coupling
    hp_ins_fsh = params["insulin_fsh_scale"] * y_ins**2 / (
        y_ins**2 + params["insulin_fsh_threshold"] ** 2
    )
    hp_p4_1 = params["hp_p4_fsh_scale"] * y_p4**2 / (
        params["hp_p4_fsh_threshold"] ** 2 + y_p4**2
    )
    hm_e2_2 = params["hm_e2_fsh_scale"] * params["hm_e2_fsh_threshold"] ** 2 / (
        params["hm_e2_fsh_threshold"] ** 2 + y_e2**2
    )
    hp_gnrh_1 = params["hp_gnrh_fsh_scale"] * y_gnrh / (
        params["hp_gnrh_fsh_threshold"] + y_gnrh
    )
    fsh_syn = params["fsh_syn_scale"] * params["hm_ih_2_threshold"] ** 5 / (
        params["hm_ih_2_threshold"] ** 5 + y_ih**5
    )
    fsh_rel = (
        params["fsh_basal_release"] + hp_p4_1 + hm_e2_2 + hp_gnrh_1
    ) * y_fshp
    f[2] = fsh_syn * hp_ins_fsh - fsh_rel
    f[3] = fsh_rel - params["fsh_clearance"] * y_fsh

    # LH with insulin coupling
    hp_ins_lh = params["insulin_lh_scale"] * y_ins**3 / (
        y_ins**3 + params["insulin_lh_threshold"] ** 3
    )
    hp_e2_2 = params["hp_e2_lh_scale"] * y_e2**2 / (
        params["hp_e2_lh_threshold"] ** 2 + y_e2**2
    )
    hm_p4_3 = params["hm_p4_lh_scale"] * params["hm_p4_lh_threshold"] ** 2 / (
        params["hm_p4_lh_threshold"] ** 2 + y_p4**2
    )
    hp_gnrh_2 = params["hp_gnrh_lh_scale"] * y_gnrh**5 / (
        params["hp_gnrh_lh_threshold"] ** 5 + y_gnrh**5
    )
    lh_syn = hp_e2_2 + hm_p4_3
    lh_rel = (params["lh_basal_release"] + hp_gnrh_2) * y_lhp
    f[4] = lh_syn * hp_ins_lh - lh_rel
    f[5] = lh_rel - params["lh_clearance"] * y_lh

    # Follicle
    hm_igf = params["igf_lh_sensitivity_scale"] * params[
        "igf_lh_sensitivity_threshold"
    ] ** 2 / (params["igf_lh_sensitivity_threshold"] ** 2 + y_igf_b**2)
    hp_lh = params["hp_lh_follicle_scale"] * y_lh**2 / (hm_igf**2 + y_lh**2)
    hp_p4_2 = params["hp_p4_follicle_scale"] * y_p4**5 / (
        params["hp_p4_follicle_threshold"] ** 5 + y_p4**5
    )
    hm_foll = params["hm_follicle_threshold"] ** 2 / (
        params["hm_follicle_threshold"] ** 2 + y_fol**2
    )
    hp_fsh_mod = params["hp_fsh_follicle_scale"] * y_fsh**2 / (
        (params["hp_fsh_follicle_threshold"] * hm_foll) ** 2 + y_fsh**2
    )
    f[6] = hp_fsh_mod - (hp_p4_2 + hp_lh) * y_fol

    # PGF, CL, P4, E2, inhibin, enzyme, OXT, IOF
    hp_enz = params["hp_enz_pg_scale"] * y_enz / (
        params["hp_enz_pg_threshold"] + y_enz
    )
    hp_ot = y_ot1**10 / (params["hp_oxt_pg_threshold"] ** 10 + y_ot1**10)
    f[7] = hp_enz * hp_ot - params["pg_clearance"] * y_pg

    hp_cl_1 = params["hp_cl_self_scale"] * y_cl**30 / (
        y_cl**30 + params["hp_cl_self_threshold"] ** 30
    )
    hp_iof = params["hp_iof_scale"] * y_iof**5 / (
        params["hp_iof_threshold"] ** 5 + y_iof**5
    )
    f[8] = params["cl_formation_scale"] * hp_lh * y_fol + hp_cl_1 - hp_iof * y_cl
    f[9] = params["p4_baseline"] + params["cl_to_p4_scale"] * y_cl**2 - params[
        "p4_clearance"
    ] * y_p4
    f[10] = params["follicle_to_e2_scale"] * y_fol**2 - params["e2_clearance"] * y_e2
    f[11] = (
        params["follicle_to_inhibin_scale"] * y_fol**2
        - params["inhibin_clearance"] * y_ih
    )
    hp_p4_3_enz = params["hp_p4_enz_scale"] * y_p4 / (
        params["hp_p4_enz_threshold"] + y_p4
    )
    f[12] = hp_p4_3_enz - params["enz_clearance"] * y_enz

    hp_e2_3 = params["hp_e2_oxt_scale"] * y_e2**2 / (
        params["hp_e2_oxt_threshold"] ** 2 + y_e2**2
    )
    input_oxt_l = (
        params["lactation_oxt_scale"] * np.exp(-params["lactation_oxt_decay"] * t**2)
        if mode == "lactating"
        else 0.0
    )
    f[13] = input_oxt_l + hp_e2_3 * y_cl**2 - params["oxt_clearance"] * y_ot1

    hp_pg = params["hp_pg_iof_scale"] * y_pg**10 / (
        y_pg**10 + params["hp_pg_iof_threshold"] ** 10
    )
    hp_cl_3 = y_cl / (y_cl + params["hp_cl_iof_threshold"])
    f[14] = hp_pg * hp_cl_3 - params["iof_clearance"] * y_iof

    # IGF
    hp_ins_igf = params["insulin_igf_scale"] * y_ins / (
        y_ins + params["insulin_igf_threshold"]
    )
    hm_p4_igf = params["hm_p4_igf_scale"] * params["hm_p4_igf_threshold"] ** 4 / (
        params["hm_p4_igf_threshold"] ** 4 + y_p4**4
    )
    f[15] = params["igf_basal_release"] + hm_p4_igf * hp_ins_igf - params[
        "igf_clearance"
    ] * y_igf_b

    # Metabolic model. Dexa multipliers are omitted because Phase 1 is no-Dexa.
    glu_pool = params["c0"] * DMI
    glu_feed_gng = (1.0 - params["feed_direct_blood_fraction"]) * glu_pool
    glu_feed_bl = params["feed_direct_blood_fraction"] * glu_pool

    ins_sec = params["insulin_secretion_max"] * y_glu**10 / (
        y_glu**10 + params["insulin_glucose_threshold"] ** 10
    )
    ins_deg = params["insulin_clearance"] * y_ins
    gluca_sec = params["glucagon_secretion_max"] * params[
        "glucagon_glucose_threshold"
    ] ** 2 / (y_glu**2 + params["glucagon_glucose_threshold"] ** 2)
    gluca_deg = params["glucagon_clearance"] * y_gluca

    glu_prod = params["hepatic_glucose_release"] * y_lv * y_gluca
    glu_bl_lv = params["blood_to_liver_max"] * y_glu**10 / (
        y_glu**10 + params["blood_to_liver_glucose_threshold"] ** 10
    ) * y_ins
    glu_lv_st = (
        params["liver_to_storage_max"]
        * (1.0 - y_glu_st / 1000.0)
        * params["milk_storage_threshold"] ** 2
        / (params["milk_storage_threshold"] ** 2 + Milk**2)
        * y_lv
        * y_ins
    )
    glu_st_lv = params["storage_to_liver_max"] * y_gluca * y_glu_st**10 / (
        y_glu_st**10 + params["storage_to_liver_threshold"] ** 10
    )
    glu_fat_lv = (
        params["fat_to_liver_max"]
        * params["fat_mobilization_storage_threshold"] ** 10
        / (y_glu_st**10 + params["fat_mobilization_storage_threshold"] ** 10)
        * y_fat
        / (y_fat + params["fat_glucose_transfer_threshold"])
        * y_gluca
    )
    glu_lv_fat = (
        params["hp_lh_follicle_threshold"]
        * params["milk_fat_threshold"]
        / (params["milk_fat_threshold"] + Milk)
        * y_glu_st**10
        / (y_glu_st**10 + params["storage_fat_threshold"] ** 10)
        * y_lv
        * y_ins
    )
    glu_bl_usage = (
        params["blood_usage_max"]
        * y_glu**10
        / (y_glu**10 + params["blood_usage_glucose_threshold"] ** 10)
        + params["milk_glucose_requirement"] * Milk
    )
    glu_lv_usage = params["liver_usage_rate"] * y_lv

    f[16] = ins_sec - ins_deg
    f[17] = (
        glu_feed_bl + glu_prod - glu_bl_lv - glu_bl_usage
    ) / params["blood_volume_l"]
    f[18] = params["fat_to_glucose"] * glu_lv_fat - glu_fat_lv
    f[19] = (
        glu_feed_gng
        - glu_prod
        + glu_st_lv
        - glu_lv_st
        + params["glucose_to_fat"] * glu_fat_lv
        - glu_lv_fat
        + glu_bl_lv
        - glu_lv_usage
    )
    f[20] = glu_lv_st - glu_st_lv
    f[21] = gluca_sec - gluca_deg

    return f
