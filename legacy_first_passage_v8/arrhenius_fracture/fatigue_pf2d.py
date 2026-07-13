"""Deprecated AT2 fatigue driver.

This project now treats fatigue crack growth through the v8 sharp-front
multifront engine, not through the AT2 phase-field relaxation driver.  The old
AT2 file could report cycle-dependent crack growth even when the Arrhenius
fatigue process-zone fields were essentially zero, because the capped AT2
iteration itself advanced the damage field.  Keeping that path runnable is too
risky.

Use one of these instead:

  python -m arrhenius_fracture.fatigue_sharp_front ...
      K-controlled V1 calibration/reduction using the same FrontEngine renewal
      and process-zone ledger.

  python -m arrhenius_fracture.sharp_front --mode 2d --crystal-aniso \
      --crystal-branch --j-decomposition cluster ...
      v8 2-D sharp-front multifront fracture/branching/deflection model.

The 2-D fatigue implementation should be added by coupling the V1 cycle-block
controller to each active v8 front.  It must not revive AT2 damage evolution as
a crack-growth criterion.
"""

from __future__ import annotations


def main(argv=None):
    raise SystemExit(
        "arrhenius_fracture.fatigue_pf2d is disabled: use the v8 sharp-front "
        "drivers (fatigue_sharp_front for K-controlled calibration, sharp_front "
        "--mode 2d for the multifront/branching model)."
    )


if __name__ == "__main__":
    main()
