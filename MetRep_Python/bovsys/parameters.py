"""MATLAB-derived parameters for the BovSys core model.

The MATLAB file ``BovSys_para_dexa_v3`` is treated as the source of truth.
Phase 1 uses parameters 1..98 for the metabolic-reproductive model and keeps
Dexa-specific PK/PD constants out of the core implementation.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Parameter:
    """One scalar model parameter with its MATLAB index and description."""

    index: int
    name: str
    value: float
    description: str


PARAMETERS: tuple[Parameter, ...] = (
    Parameter(1, "gnrh_hypo_max", 16.0, "maximum hypothalamic GnRH store"),
    Parameter(2, "gnrh_syn", 2.75, "GnRH synthesis rate"),
    Parameter(3, "gnrh_release_scale", 2.05, "GnRH release scale"),
    Parameter(4, "hm_e2_1_threshold", 0.0972, "E2 inhibition threshold for GnRH"),
    Parameter(5, "hm_p4_1_threshold", 0.35, "P4 inhibition threshold for GnRH"),
    Parameter(6, "hm_p4_2_scale", 1.91, "P4-related GnRH release component"),
    Parameter(7, "hm_p4_2_threshold", 0.252, "P4 threshold for GnRH release component"),
    Parameter(8, "hp_e2_1_scale", 0.99, "E2 positive GnRH gate scale"),
    Parameter(9, "hp_e2_1_threshold", 0.648, "E2 positive GnRH gate threshold"),
    Parameter(10, "gnrh_clearance", 1.63, "GnRH clearance"),
    Parameter(11, "fsh_syn_scale", 4.21, "FSH synthesis scale"),
    Parameter(12, "hm_ih_2_threshold", 0.118, "inhibin threshold for FSH synthesis"),
    Parameter(13, "hp_p4_fsh_scale", 0.293, "P4 effect on FSH release"),
    Parameter(14, "hp_p4_fsh_threshold", 0.152, "P4 threshold for FSH release"),
    Parameter(15, "hm_e2_fsh_scale", 0.396, "E2 effect on FSH release"),
    Parameter(16, "hm_e2_fsh_threshold", 0.312, "E2 threshold for FSH release"),
    Parameter(17, "hp_gnrh_fsh_scale", 1.23, "GnRH effect on FSH release"),
    Parameter(18, "hp_gnrh_fsh_threshold", 0.0708, "GnRH threshold for FSH release"),
    Parameter(19, "fsh_clearance", 2.73, "blood FSH clearance"),
    Parameter(20, "fsh_basal_release", 0.948, "basal FSH release"),
    Parameter(21, "hp_e2_lh_scale", 0.376, "E2 effect on LH synthesis"),
    Parameter(22, "hp_e2_lh_threshold", 0.243, "E2 threshold for LH synthesis"),
    Parameter(23, "hm_p4_lh_scale", 2.71, "P4 effect on LH synthesis"),
    Parameter(24, "hm_p4_lh_threshold", 0.0269, "P4 threshold for LH synthesis"),
    Parameter(25, "hp_gnrh_lh_scale", 2.22, "GnRH effect on LH release"),
    Parameter(26, "hp_gnrh_lh_threshold", 0.69, "GnRH threshold for LH release"),
    Parameter(27, "lh_basal_release", 0.0141, "basal LH release"),
    Parameter(28, "lh_clearance", 2.0, "blood LH clearance"),
    Parameter(29, "hp_fsh_follicle_scale", 0.562, "FSH-driven follicle growth scale"),
    Parameter(30, "hp_fsh_follicle_threshold", 1.497, "FSH threshold for follicle growth"),
    Parameter(31, "hm_follicle_threshold", 0.322, "follicle self-limitation threshold"),
    Parameter(32, "hp_p4_follicle_scale", 1.1, "P4-driven follicle loss scale"),
    Parameter(33, "hp_p4_follicle_threshold", 0.126, "P4 threshold for follicle loss"),
    Parameter(34, "hp_lh_follicle_scale", 3.49, "LH-driven follicle loss/ovulation scale"),
    Parameter(35, "hp_lh_follicle_threshold", 0.1, "also used by MATLAB metabolic liver-to-fat flux"),
    Parameter(36, "hp_enz_pg_scale", 53.91, "enzyme effect on PGF production"),
    Parameter(37, "hp_enz_pg_threshold", 1.43, "enzyme threshold for PGF production"),
    Parameter(38, "hp_oxt_pg_threshold", 1.087, "OXT threshold for PGF production"),
    Parameter(39, "pg_clearance", 1.23, "PGF clearance"),
    Parameter(40, "cl_formation_scale", 0.4, "follicle-to-CL formation scale"),
    Parameter(41, "hp_cl_self_scale", 0.0335, "CL self-support scale"),
    Parameter(42, "hp_cl_self_threshold", 0.2807, "CL self-support threshold"),
    Parameter(43, "hp_iof_scale", 41.39, "IOF luteolysis scale"),
    Parameter(44, "hp_iof_threshold", 1.82, "IOF luteolysis threshold"),
    Parameter(45, "cl_to_p4_scale", 2.25, "CL-dependent P4 production"),
    Parameter(46, "p4_clearance", 1.41, "P4 clearance"),
    Parameter(47, "p4_baseline", 0.1, "basal P4 production"),
    Parameter(48, "follicle_to_e2_scale", 2.19, "follicle-dependent E2 production"),
    Parameter(49, "e2_clearance", 1.23, "E2 clearance"),
    Parameter(50, "follicle_to_inhibin_scale", 1.41, "follicle-dependent inhibin production"),
    Parameter(51, "inhibin_clearance", 0.475, "inhibin clearance"),
    Parameter(52, "hp_p4_enz_scale", 3.58, "P4 effect on enzyme production"),
    Parameter(53, "hp_p4_enz_threshold", 0.77, "P4 threshold for enzyme production"),
    Parameter(54, "enz_clearance", 2.98, "enzyme clearance"),
    Parameter(55, "hp_e2_oxt_scale", 1.59, "E2/CL effect on OXT production"),
    Parameter(56, "hp_e2_oxt_threshold", 0.143, "E2 threshold for OXT production"),
    Parameter(57, "oxt_clearance", 0.644, "OXT clearance"),
    Parameter(58, "lactation_oxt_scale", 1.5, "lactation OXT input scale"),
    Parameter(59, "fat_mobilization_storage_threshold", 10.0, "storage threshold in fat-to-liver flux"),
    Parameter(60, "lactation_oxt_decay", 0.0007, "lactation OXT input decay"),
    Parameter(61, "hp_pg_iof_scale", 39.68, "PGF effect on IOF production"),
    Parameter(62, "hp_pg_iof_threshold", 1.22, "PGF threshold for IOF production"),
    Parameter(63, "hp_cl_iof_threshold", 0.6, "CL gate threshold for IOF production"),
    Parameter(64, "iof_clearance", 0.298, "IOF clearance"),
    Parameter(65, "hm_p4_igf_scale", 2.5, "P4 modulation of IGF drive"),
    Parameter(66, "hm_p4_igf_threshold", 0.3, "P4 threshold for IGF drive"),
    Parameter(67, "igf_clearance", 1.7, "IGF clearance"),
    Parameter(68, "igf_basal_release", 0.4, "basal IGF release"),
    Parameter(69, "igf_lh_sensitivity_scale", 1.0, "IGF modulation of LH follicle sensitivity"),
    Parameter(70, "igf_lh_sensitivity_threshold", 0.5, "IGF threshold for LH follicle sensitivity"),
    Parameter(71, "insulin_fsh_scale", 3.0, "insulin effect on FSH"),
    Parameter(72, "insulin_fsh_threshold", 15.0, "insulin threshold for FSH"),
    Parameter(73, "insulin_lh_scale", 1.05, "insulin effect on LH"),
    Parameter(74, "insulin_lh_threshold", 16.0, "insulin threshold for LH"),
    Parameter(75, "insulin_igf_scale", 0.4, "insulin effect on IGF"),
    Parameter(76, "insulin_igf_threshold", 15.0, "insulin threshold for IGF"),
    Parameter(77, "feed_direct_blood_fraction", 0.08, "fraction of feed glucose routed directly to blood"),
    Parameter(78, "fat_to_liver_max", 8.0, "fat to liver glucose-equivalent transfer maximum"),
    Parameter(79, "storage_to_liver_threshold", 10.0, "storage threshold for storage-to-liver flux"),
    Parameter(80, "storage_to_liver_max", 1350.0, "storage to liver transfer maximum"),
    Parameter(81, "milk_storage_threshold", 10.0, "milk threshold for liver-to-storage flux"),
    Parameter(82, "liver_to_storage_max", 180.0, "liver to storage transfer maximum"),
    Parameter(83, "blood_to_liver_glucose_threshold", 0.45, "blood glucose threshold for blood-to-liver flux"),
    Parameter(84, "glucagon_secretion_max", 70182.0, "glucagon secretion scale"),
    Parameter(85, "glucagon_clearance", 350.87, "glucagon clearance"),
    Parameter(86, "glucagon_glucose_threshold", 0.5, "glucose threshold for glucagon secretion"),
    Parameter(87, "hepatic_glucose_release", 0.0684, "glucagon-driven liver glucose release scale"),
    Parameter(88, "blood_to_liver_max", 50.0, "blood-to-liver insulin-dependent transfer scale"),
    Parameter(89, "milk_fat_threshold", 10.0, "milk threshold for liver-to-fat flux"),
    Parameter(90, "storage_fat_threshold", 1000.0, "storage threshold for liver-to-fat flux"),
    Parameter(91, "insulin_glucose_threshold", 0.5, "glucose threshold for insulin secretion"),
    Parameter(92, "blood_usage_glucose_threshold", 0.5, "glucose threshold for non-mammary usage"),
    Parameter(93, "blood_usage_max", 1000.0, "non-mammary blood glucose usage scale"),
    Parameter(94, "liver_usage_rate", 5.0, "liver glucose usage rate"),
    Parameter(95, "milk_glucose_requirement", 72.0, "glucose required per unit milk"),
    Parameter(96, "insulin_secretion_max", 84211.0, "insulin secretion scale"),
    Parameter(97, "insulin_clearance", 2105.0, "insulin clearance"),
    Parameter(98, "fat_glucose_transfer_threshold", 150000.0, "fat saturation threshold"),
)


def default_parameters() -> dict[str, float]:
    """Return a mutable name-to-value parameter dictionary."""

    params = {parameter.name: parameter.value for parameter in PARAMETERS}
    params["c0"] = 0.08
    params["fat_to_glucose"] = 9.3 / 4.1
    params["glucose_to_fat"] = 4.1 / 9.3
    params["blood_volume_l"] = 22.8
    return params


def matlab_parameter_vector() -> list[float]:
    """Return the 1..98 MATLAB parameter vector as a Python list."""

    return [parameter.value for parameter in PARAMETERS]


def parameter_table() -> list[dict[str, float | int | str]]:
    """Return parameter metadata suitable for saving as a table."""

    return [
        {
            "matlab_index": parameter.index,
            "name": parameter.name,
            "value": parameter.value,
            "description": parameter.description,
        }
        for parameter in PARAMETERS
    ]
