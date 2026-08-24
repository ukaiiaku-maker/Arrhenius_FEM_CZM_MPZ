# One D Driving Force Semantics Audit

Autonomous v9.13 computes `K_i(U)=K_per_U[i] U`. The map was derived from pre-event `KJ/U` values of an earlier 2-D reference trajectory, so it is an event-indexed reduced-mechanics map. It is neither a compliance-derived finite-specimen G nor an applied remote K. Native PF 1-D instead prescribes `Kdot*time`; it also has no structural compliance-derived G.

V2 separates the event-induced change `K_{i+1}(U_i)-K_i(U_i)` from reload `K_{i+1}(U_{i+1})-K_{i+1}(U_i)`. No pre-event field is compared to post-event geometry without naming the phase. Legacy event-wise endpoint/rise columns remain available, explicitly deprecated as resistance measures.
