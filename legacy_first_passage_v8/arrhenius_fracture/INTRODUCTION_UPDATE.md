# Introduction update — sharp-front realization and the regime span

Drop-in revisions to "Brittle Fracture as an Arrhenius Instability". The existing
intro already seeds shielding via sigma -> sigma - sigma_back and frames fracture
as a first-event competition of additive hazards. Two things are missing and
should be added: (i) a clear statement of the *specific* competing channels and
their variables in the crack-advance hazard, and (ii) the new physical result
that the *temperature dependence* of toughness — DBTT vs weakly-T vs
ceramic-softening — is selected by the sign and saturation of the
plasticity->cleavage coupling, not by an imposed switch.

## Variable table (define once, near the Arrhenius-fracture formulation)

| symbol | meaning | role |
|--------|---------|------|
| K, sigma_tip = K/sqrt(2 pi r_eff) | nominal SIF, near-tip opening stress at process-zone radius r_eff | drives the cleavage hazard |
| Lambda = ∫ lambda dt | cumulative hazard; fracture onset at Lambda ~ 1 | first-passage criterion |
| lambda_cleave = nu0 exp(-Geff/kT) | crack-advance hazard (multihit first passage over tau_c, m bonds) | the event being raced |
| Geff = H(sigma_eff) - T S - dG_emb | effective barrier of crack advance | the contested quantity |
| **sigma_eff = sigma_tip - chi_shield · sigma_back** | shielded driving stress | **toughening channel** (RAISES Geff; bounded above at H0 where sigma_eff->0) |
| sigma_back ∝ N_em | crack-tip pile-up back stress | shielding magnitude |
| **dG_emb ∝ N_em** | stored-energy / chemical barrier reduction | **embrittling channel** (LOWERS Geff) |
| N_em | emitted-dislocation ledger (grows with T as emission activates) | common driver of BOTH channels |
| **N_sat, k_rec** | saturation density / dynamic-recovery rate of N_em | **bounds both channels (storage vs recovery)** |
| H0 | zero-stress cleavage barrier | shelf height / ceiling on shielding |
| S (negative) | activation entropy of cleavage | intrinsic softening with T |

chi_shield is the dimensionless scalar form of the introduction's sigma ->
sigma - sigma_back substitution applied to the crack-opening hazard. N_sat / k_rec
encode the dislocation storage-recovery balance (Kocks-Mecking-type), absent from
the minimal model and required for a physically bounded back stress and stored
energy.

## New passage 1 — the two competing channels (after the sigma->sigma-sigma_back paragraph)

> Crack-tip plasticity couples to the crack-advance hazard through two channels of
> opposite sign, both carried by the same population of emitted dislocations N_em.
> The pile-up back stress sigma_back enters as a *shielding* term, sigma_eff =
> sigma_tip - chi_shield sigma_back, which raises the activation barrier and is the
> only channel that strengthens the tip as plasticity develops. The stored elastic
> energy (and any chemically assisted bond weakening it enables) enters as an
> *embrittling* term dG_emb that lowers the barrier. Because dislocation emission is
> itself thermally activated, both channels switch on with temperature; the
> temperature dependence of the apparent toughness is therefore set by which channel
> dominates, and by whether that dominance persists as T increases.

## New passage 2 — regime selection (replaces any single-mechanism DBTT statement)

> Three qualitatively different temperature dependences emerge from this single
> competition without any imposed transition rule. When embrittlement dominates, the
> barrier falls with temperature and the toughness decreases monotonically — the
> thermally-activated subcritical-growth behavior characteristic of oxide ceramics
> and glasses. When shielding dominates and remains effective, the barrier is held
> up as plasticity develops and the toughness rises with temperature — a ductile-to-
> brittle transition. The balance between the two gives a weakly temperature-
> dependent toughness. A fourth, non-monotonic response — a toughness peak — appears
> when shielding wins at intermediate temperature but embrittlement reclaims the tip
> at higher temperature. The selecting variables are (i) the shielding strength
> chi_shield, which sets the onset/low-temperature side of the transition, and (ii)
> the survival of a tough upper shelf, set by the saturated stored energy relative
> to the bond barrier.

## New passage 3 — the saturation/recovery requirement (a genuine model prediction)

> A specific and falsifiable consequence follows. Because both shielding and
> embrittlement scale with the dislocation ledger, an *unbounded* ledger guarantees
> that embrittlement — which lowers the barrier without bound — wins the high-
> temperature limit, so no tough upper shelf can exist. A ductile-to-brittle
> transition with a stable upper shelf therefore *requires* the stored-energy
> embrittlement to saturate, which in turn requires the underlying dislocation
> density to saturate through dynamic recovery (storage balanced by annihilation).
> In the model this enters as a recovery term on N_em; the saturated value of dG_emb
> relative to H0 then determines whether the high-temperature limit is tough (DBTT)
> or brittle (ceramic). The recovery rate, which is independently constrained by the
> plasticity law, thus becomes a predictor of upper-shelf toughness.

## Caveats to state honestly in the text (revised after testing)

> The model produces a finite, two-shelf ductile-to-brittle transition: a lower
> (brittle) shelf, a graded rise spanning a few tens of kelvin, and a finite upper
> (tough) shelf whose height is set by the stress at which the rising near-tip field
> overcomes the saturated shielding back stress. (An earlier apparent "arrest" upper
> limit was a measurement artifact of a loading window set below the upper shelf, not
> a feature of the model.) The post-initiation crack velocity saturates at the
> Rayleigh speed through a finite v(K) law, recovering the Wiederhorn Region-I/III
> picture; this bounds unstable fast fracture but does not affect the first-passage
> toughness Kc. The width and position of the transition are governed by the
> temperature dependence of the emission entropy (the activation-entropy crossover);
> a smoother entropy crossover broadens the transition and the associated toughness
> peak. These are refinements of the kinetic detail, not of the regime-selection
> mechanism, which is fixed by the shielding/embrittlement competition and the
> saturation of the dislocation ledger.
