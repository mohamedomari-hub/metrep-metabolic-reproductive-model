"""Feeding and lactation scenario definitions."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
METREP_FORCING_CSV = PROJECT_ROOT / "data" / "DM_MILK_Data.csv"

SCENARIO_DEFAULT_DAYS = {
    "baseline_non_lactating": 60.0,
    "acute_negative_energy_balance": 90.0,
    "non_lactating_acute": 90.0,
    "chronic_negative_energy_balance": 330.0,
    "non_lactating_chronic": 330.0,
    "lactating_c0_20": 271.1891,
    "lactating_c0_22_5": 271.1891,
    "lactating_c0_25": 271.1891,
    "lactating_c0_30": 271.1891,
}


@dataclass(frozen=True)
class Scenario:
    """Exogenous forcing definition for a simulation."""

    name: str
    mode: str
    c0: float
    t_eval: np.ndarray
    dmi: np.ndarray
    milk: np.ndarray
    description: str


def time_grid(days: float, step: float = 1.0 / 48.0) -> np.ndarray:
    """Return an inclusive simulation grid in days."""

    days = float(days)
    t = np.arange(0.0, days, step, dtype=float)
    if t.size == 0 or not np.isclose(t[-1], days):
        t = np.append(t, days)
    return t


def constant_non_lactating(days: float = 60.0, step: float = 1.0 / 48.0) -> Scenario:
    """MATLAB non-lactating standard feeding: 11.7 kg DMI/day, no milk."""

    t = time_grid(days, step)
    return Scenario(
        name="non_lactating_standard",
        mode="non_lactating",
        c0=0.08,
        t_eval=t,
        dmi=np.full_like(t, 11700.0),
        milk=np.zeros_like(t),
        description="Non-lactating baseline: constant DMI 11700 g/day, milk 0.",
    )


def constant_lactating(
    days: float = 60.0,
    step: float = 1.0 / 48.0,
    dmi_g_per_day: float = 20000.0,
    milk_l_per_day: float = 30.0,
    c0: float = 0.20,
) -> Scenario:
    """Simple lactating scenario with constant DMI and milk yield."""

    t = time_grid(days, step)
    return Scenario(
        name="constant_lactating",
        mode="lactating",
        c0=c0,
        t_eval=t,
        dmi=np.full_like(t, dmi_g_per_day),
        milk=np.full_like(t, milk_l_per_day),
        description=(
            "Lactating constant forcing: "
            f"DMI {dmi_g_per_day:g} g/day, milk {milk_l_per_day:g} L/day."
        ),
    )


def early_lactation_curve(
    days: float = 120.0,
    step: float = 1.0 / 48.0,
    peak_milk_l_per_day: float = 35.0,
    peak_day: float = 45.0,
    start_milk_l_per_day: float = 18.0,
    end_milk_l_per_day: float = 24.0,
    start_dmi_g_per_day: float = 14000.0,
    peak_dmi_g_per_day: float = 22000.0,
    c0: float = 0.20,
) -> Scenario:
    """Approximate early-lactation rise and decline for DMI and milk.

    This is a practical built-in forcing curve, not a Dexa scenario. Use
    :func:`scenario_from_csv` when exact paper or PDF forcing values are known.
    """

    t = time_grid(days, step)
    peak_day = min(max(float(peak_day), 0.0), float(days))
    milk = np.interp(
        t,
        [0.0, peak_day, float(days)],
        [start_milk_l_per_day, peak_milk_l_per_day, end_milk_l_per_day],
    )
    dmi = np.interp(
        t,
        [0.0, peak_day, float(days)],
        [start_dmi_g_per_day, peak_dmi_g_per_day, 0.95 * peak_dmi_g_per_day],
    )
    return Scenario(
        name="early_lactation_curve",
        mode="lactating",
        c0=c0,
        t_eval=t,
        dmi=dmi,
        milk=milk,
        description=(
            "Approximate early lactation forcing: milk rises to "
            f"{peak_milk_l_per_day:g} L/day at day {peak_day:g}; DMI rises to "
            f"{peak_dmi_g_per_day:g} g/day."
        ),
    )


def negative_energy_balance(days: float = 90.0, step: float = 1.0 / 48.0) -> Scenario:
    """Acute restriction from MATLAB: DMI reduced to 33% for 15 days."""

    t = time_grid(days, step)
    dmi = np.full_like(t, 11700.0)
    mask = (t >= 46.0) & (t <= 61.0)
    dmi[mask] = 0.33 * 11700.0
    return Scenario(
        name="non_lactating_acute",
        mode="non_lactating",
        c0=0.08,
        t_eval=t,
        dmi=dmi,
        milk=np.zeros_like(t),
        description="Acute negative energy balance: 33% DMI from day 46 to 61.",
    )


def chronic_negative_energy_balance(days: float = 330.0, step: float = 1.0 / 48.0) -> Scenario:
    """MetRep chronic restriction/refeeding: 58% DMI for 210 days, then 70 days refeeding."""

    t = time_grid(days, step)
    dmi = np.full_like(t, 11700.0)
    restriction = (t >= 45.0) & (t <= 45.0 + 210.0)
    refeeding = (t >= 45.0 + 210.0) & (t <= 45.0 + 210.0 + 70.0)
    dmi[restriction] = 0.58 * 11700.0
    dmi[refeeding] = 1.6 * 11700.0
    return Scenario(
        name="non_lactating_chronic",
        mode="non_lactating",
        c0=0.08,
        t_eval=t,
        dmi=dmi,
        milk=np.zeros_like(t),
        description=(
            "MetRep chronic negative energy balance: 58% DMI from day 45 "
            "for 210 days, then 160% DMI refeeding for 70 days."
        ),
    )


def improved_feeding(days: float = 90.0, step: float = 1.0 / 48.0) -> Scenario:
    """Improved feeding scenario for comparison with restriction."""

    t = time_grid(days, step)
    dmi = np.full_like(t, 11700.0)
    dmi[t >= 46.0] = 1.25 * 11700.0
    return Scenario(
        name="improved_feeding",
        mode="non_lactating",
        c0=0.08,
        t_eval=t,
        dmi=dmi,
        milk=np.zeros_like(t),
        description="Improved feeding: DMI increased by 25% from day 46.",
    )


def lactating_metrep_forcing(
    days: float | None = None,
    step: float = 1.0 / 48.0,
    c0: float = 0.20,
    name: str = "lactating_c0_20",
) -> Scenario:
    """MetRep lactating scenario from weekly DMI/milk forcing data."""

    table = pd.read_csv(METREP_FORCING_CSV, header=None, names=["DMI", "Milk", "weeks"])
    source_time = table["weeks"].to_numpy(dtype=float) * 7.0
    if days is None:
        days = float(source_time[-1])
    t = time_grid(days, step)
    dmi = np.interp(t, source_time, table["DMI"].to_numpy(dtype=float))
    milk = np.interp(t, source_time, table["Milk"].to_numpy(dtype=float))
    return Scenario(
        name=name,
        mode="lactating",
        c0=c0,
        t_eval=t,
        dmi=dmi,
        milk=milk,
        description=(
            f"MetRep lactating scenario from weekly DMI/milk data, c0={c0:g}; "
            f"default horizon {source_time[-1]:.4g} days."
        ),
    )


def lactating_c0_20(days: float | None = None, step: float = 1.0 / 48.0) -> Scenario:
    """MetRep lactating scenario with c0=20%."""

    return lactating_metrep_forcing(days=days, step=step, c0=0.20, name="lactating_c0_20")


def lactating_c0_22_5(days: float | None = None, step: float = 1.0 / 48.0) -> Scenario:
    """MetRep lactating scenario with c0=22.5%."""

    return lactating_metrep_forcing(days=days, step=step, c0=0.225, name="lactating_c0_22_5")


def lactating_c0_25(days: float | None = None, step: float = 1.0 / 48.0) -> Scenario:
    """MetRep lactating scenario with c0=25%."""

    return lactating_metrep_forcing(days=days, step=step, c0=0.25, name="lactating_c0_25")


def lactating_c0_30(days: float | None = None, step: float = 1.0 / 48.0) -> Scenario:
    """MetRep lactating scenario with c0=30%."""

    return lactating_metrep_forcing(days=days, step=step, c0=0.30, name="lactating_c0_30")


def custom_step_feeding(
    days: float = 90.0,
    step: float = 1.0 / 48.0,
    mode: str = "non_lactating",
    baseline_dmi_g_per_day: float = 11700.0,
    baseline_milk_l_per_day: float = 0.0,
    dmi_multiplier: float = 1.0,
    start_day: float = 46.0,
    end_day: float | None = None,
    c0: float | None = None,
    name: str = "custom_step_feeding",
) -> Scenario:
    """Build a scenario with a step change in DMI over a time interval."""

    t = time_grid(days, step)
    dmi = np.full_like(t, baseline_dmi_g_per_day)
    milk = np.full_like(t, baseline_milk_l_per_day)
    if end_day is None:
        mask = t >= start_day
    else:
        mask = (t >= start_day) & (t <= end_day)
    dmi[mask] = dmi_multiplier * baseline_dmi_g_per_day
    if c0 is None:
        c0 = 0.20 if mode == "lactating" else 0.08
    return Scenario(
        name=name,
        mode=mode,
        c0=float(c0),
        t_eval=t,
        dmi=dmi,
        milk=milk,
        description=(
            f"{mode} step feeding: DMI x{dmi_multiplier:g} from day "
            f"{start_day:g}" + ("" if end_day is None else f" to {end_day:g}") + "."
        ),
    )


def lactating_from_series(
    t_eval: np.ndarray,
    dmi: np.ndarray,
    milk: np.ndarray,
    c0: float = 0.20,
    name: str = "lactating",
) -> Scenario:
    """Build a lactating scenario from externally loaded DMI and milk vectors."""

    return Scenario(
        name=name,
        mode="lactating",
        c0=c0,
        t_eval=np.asarray(t_eval, dtype=float),
        dmi=np.asarray(dmi, dtype=float),
        milk=np.asarray(milk, dtype=float),
        description="Lactating scenario from external forcing data.",
    )


def scenario_from_csv(
    csv_path: str,
    step: float = 1.0 / 48.0,
    mode: str = "lactating",
    c0: float | None = None,
    name: str | None = None,
    time_column: str = "time_days",
    dmi_column: str = "DMI",
    milk_column: str = "Milk",
) -> Scenario:
    """Build a scenario from CSV columns for time, DMI, and milk.

    The CSV may use sparse time points; DMI and milk are linearly interpolated
    onto a regular grid for the solver.
    """

    table = pd.read_csv(csv_path)
    missing = [col for col in [time_column, dmi_column, milk_column] if col not in table.columns]
    if missing:
        raise ValueError(f"Missing required CSV columns: {', '.join(missing)}")

    source_time = table[time_column].to_numpy(dtype=float)
    order = np.argsort(source_time)
    source_time = source_time[order]
    if source_time.size < 2:
        raise ValueError("Scenario CSV must contain at least two time points.")
    if np.any(np.diff(source_time) <= 0):
        raise ValueError("Scenario CSV time values must be unique and increasing.")

    t = time_grid(float(source_time[-1]), step)
    dmi = np.interp(t, source_time, table[dmi_column].to_numpy(dtype=float)[order])
    milk = np.interp(t, source_time, table[milk_column].to_numpy(dtype=float)[order])
    if c0 is None:
        c0 = 0.20 if mode == "lactating" else 0.08

    return Scenario(
        name=name or str(csv_path).split("/")[-1].rsplit(".", 1)[0],
        mode=mode,
        c0=float(c0),
        t_eval=t,
        dmi=dmi,
        milk=milk,
        description=f"{mode} scenario loaded from {csv_path}.",
    )


def feeding_scenarios(days: float | None = None, step: float = 1.0 / 48.0) -> list[Scenario]:
    """Return MetRep non-lactating feeding scenarios."""

    return [
        built_in_scenario("baseline_non_lactating", days=days, step=step),
        built_in_scenario("acute_negative_energy_balance", days=days, step=step),
        built_in_scenario("chronic_negative_energy_balance", days=days, step=step),
    ]


SCENARIO_BUILDERS = {
    "baseline_non_lactating": constant_non_lactating,
    "acute_negative_energy_balance": negative_energy_balance,
    "chronic_negative_energy_balance": chronic_negative_energy_balance,
    "lactating_c0_20": lactating_c0_20,
    "lactating_c0_22_5": lactating_c0_22_5,
    "lactating_c0_25": lactating_c0_25,
    "lactating_c0_30": lactating_c0_30,
}

SCENARIO_ALIASES = {
    "negative_energy_balance": "acute_negative_energy_balance",
    "non_lactating_acute": "acute_negative_energy_balance",
    "non_lactating_chronic": "chronic_negative_energy_balance",
}


def available_scenarios() -> list[str]:
    """Return names of built-in scenario constructors."""

    return list(SCENARIO_BUILDERS)


def scenario_default_days(name: str) -> float:
    """Return the MetRep default horizon for a built-in scenario."""

    name = SCENARIO_ALIASES.get(name, name)
    if name not in SCENARIO_DEFAULT_DAYS:
        raise ValueError(f"No default horizon defined for scenario {name!r}.")
    return SCENARIO_DEFAULT_DAYS[name]


def built_in_scenario(name: str, days: float | None = None, step: float = 1.0 / 48.0) -> Scenario:
    """Return a built-in scenario by name."""

    name = SCENARIO_ALIASES.get(name, name)
    if name not in SCENARIO_BUILDERS:
        raise ValueError(f"Unknown scenario {name!r}. Available: {', '.join(available_scenarios())}")
    if days is None:
        days = scenario_default_days(name)
    return SCENARIO_BUILDERS[name](days=days, step=step)
