"""Dexa scenario extension placeholder.

Phase 1 deliberately keeps Dexa outside the core metabolic-reproductive model.
The MATLAB Dexa-specific components are:
    - states 23..25: A_dep, A_cent, C_e
    - PK constants: ka, ke, F, keo
    - PD multipliers: Effect_gluca and Effect_bt
    - dosing logic in solve_with_dose
"""


def main() -> None:
    raise NotImplementedError(
        "Dexa scenarios are intentionally deferred until the 22-state "
        "metabolic-reproductive core is validated against MATLAB outputs."
    )


if __name__ == "__main__":
    main()
