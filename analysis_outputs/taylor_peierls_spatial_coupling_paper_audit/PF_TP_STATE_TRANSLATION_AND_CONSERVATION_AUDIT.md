# PF Taylor/Peierls state translation and conservation audit

The moving process zone is translated continuously inside the Strang-coupled kinetic integrator; the outer `n_fire` commit does not apply a second population reset. The neutral observer serializes accepted-row endpoints, not every inner microstep. Therefore `pf_tp_event_state_translation.csv` distinguishes the exact accepted-row population balance from an exact application of the production moving-frame operator to the prior accepted profile.

For every one of the 472 physical transactions, the operator-projected conservation residual closes to at most 2.328e-10 line-count units. The accepted-row source ledger closes to 0.000e+00 after explicitly accounting for emission, escape, recovery, and inferred finite-wake discard/unobserved removal.

Summed over the 59 transactions in each candidate, the exact moving-frame operator deposits 8.08e-08–5.34e-07 mobile and 1.66e+06–1.66e+06 retained line-count units into the wake. Finite-window discard totals 3.75e+06–3.75e+06. The accepted-row brackets add 7.3e+05–7.41e+05 by emission, report 1.29e-06–3.37e-06 escaped and 0–0 recovered, and renew no persistent source sites. These sums are flow-through totals, not the final stored population.

At the two physical onset boundaries, regional columns report the active 0–2 um and 0–10 um state and wake bands 0–10, 10–50, 50–100, and beyond 100 um. The production wake support ends at 100 um, so the farther-wake row is identically zero; older state is recorded as discarded rather than silently retained. Candidate diversity at re-initiation resides chiefly in mobile-versus-retained partitioning over the active/wake ledger. The exponentially weighted near-tip unsigned sum is nearly invariant.

No stochastic trajectory was replayed. `PREVIOUS_ACCEPTED_ROW_TO_EVENT_CROSSING_ROW` is the strongest time bracket supported by the immutable observer; it is not mislabelled as an inner-microstep snapshot.
