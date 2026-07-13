# R-curve-like outputs for the four-class 500 um sweep

Each completed class-temperature case retains the full raw solver history in
`steps_<T>K.csv` and additionally writes:

- `R_curve_event_sampled.csv`: accepted crack-growth events only;
- `R_curve_K_vs_crack_extension.png`: KJ versus projected crack extension.

The event-sampled CSV contains, when available:

- growth_event_id
- step
- crack_extension_um
- a_tip_mm
- da_block_um
- KJ_MPa_sqrt_m
- N_em
- B
- n_fire
- sigma_tip_GPa
- sigma_back_GPa
- lambda_c_per_s
- lambda_e_per_s
- G_cleave_eff_eV

The top-level `four_class_temperature_summary.csv` also reports late-growth
propagation metrics over 200 <= delta-a_x <= 500 um:

- Kprop_200_500um_median
- Kprop_200_500um_mean
- Kprop_200_500um_p10
- Kprop_200_500um_p90
- delta_KR_median_minus_init
- n_growth_events

Class directories contain `R_curves_all_temperatures.png`, and the sweep root
contains `four_class_init_vs_propagation_K_vs_T.png`.

These are described as R-curve-like propagation-resistance outputs because the
x-axis is projected crack extension and the ordinate is the front-specific KJ
sampled at accepted crack-growth events; they are not asserted to be a
standardized ASTM R-curve measurement.
