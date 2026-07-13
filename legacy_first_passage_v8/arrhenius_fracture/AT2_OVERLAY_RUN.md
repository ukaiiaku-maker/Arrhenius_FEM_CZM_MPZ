# AT2 overlay -- how to run

Run with `-m` from the directory that CONTAINS the `arrhenius_fracture/` folder.
Each command below is a SINGLE line -- copy the line only, not the description.

## Temperature sweep
python3 -m arrhenius_fracture.at2_overlay --mode sweep --temperatures 500 650 800 950 1100 --out sweep_out

## Single demo
python3 -m arrhenius_fracture.at2_overlay --mode demo --out demo_out

## Propagation-only baseline (clean, fully reproducible)
python3 -m arrhenius_fracture.at2_overlay --mode demo --nucleation-off --out demo_nonucl

## Lower load to find arrest / DBTT-like behavior
python3 -m arrhenius_fracture.at2_overlay --mode sweep --temperatures 500 650 800 950 1100 --sigma0 0.15e9 --dsigma 0.015e9 --out sweep_lowload

Outputs:
  sweep_out/at2_sweep.png        one phi panel per T + response-vs-T plot
  sweep_out/sweep_summary.json   main_tip_mm, final_cracks, total_nucleations, damaged_area_mm2, resid_median
  demo_out/at2_overlay.png       phi evolution, Gc microstructure, advance handshake, crack inventory

Flags: --steps --sigma0 --dsigma --dt --nx --ny --Lx --ell --seed --nucleation-off --nu0 --dG0-eV --sigma-star-GPa --no-render

## Notes
- main crack now MEANDERS (follows low-Gc paths) and propagation is REPRODUCIBLE
  across machines (tip direction no longer depends on the skimage skeletonize
  version). Run with --nucleation-off for the cleanest single-crack baseline.
- nucleation PLACEMENT and DENSITY are still only qualitative: the reduced stress
  proxy lacks the true spatial stress structure, so secondary cracks seed at the
  weakest microstructure spots (a boundary band keeps them off the very edges, but
  they are not yet correctly concentrated ahead of the tip). This is the [FEM] item
  -- the real stress field fixes placement. Read sweeps for TRENDS vs T/seed, not
  absolute counts.
