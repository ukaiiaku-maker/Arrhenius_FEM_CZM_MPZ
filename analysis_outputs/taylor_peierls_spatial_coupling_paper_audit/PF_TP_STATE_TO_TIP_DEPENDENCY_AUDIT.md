# PF Taylor/Peierls state-to-tip dependency audit

## Production source conclusion

The eight cases use the v10.2.21 persistent-site source with the v10.2.22 mesh-independent physical along-front-width overlay, over the v10.2.14 active-only signed shielding atlas. The production fracture law does **not** consume total or retained-wake content directly. It consumes a near-tip projection of the active state through radius, density-limited source width/multiplicity, unsigned Taylor backstress, and active signed shielding. The retained wake remains serialized and advected, but its measured shielding operator is identically zero and it is absent from backstress, radius, and multiplicity.

## Direct answers

1. Backstress uses all 80 active bins with `exp[-x/max(0.5 um, dx)]` weighting of unsigned mobile plus retained content. Here `dx=0.625 um`, so the actual weighting length is 0.625 um and the first bins dominate.
2. Shielding uses signed retained content in all 80 active bins and the extension-dependent measured physical-x kernel. Mobile shielding is zero.
3. Retained material behind the advanced tip remains in the wake ledger but is mechanically decoupled in this stack: wake shielding is disabled and wake state does not enter local source geometry.
4. Backstress and shielding are evaluated in tip-relative coordinates. The tensor probe begins in the laboratory FEM mesh and is projected into the current tip frame.
5. Active mobile, retained, signed species, and accumulated slip translate toward the moving tip.
6. Fractions crossing the tip are deposited in the laboratory-frame wake; old wake advects farther behind the new tip.
7. Persistent source sites are neither depleted nor refreshed. Population beyond the 100-um wake support is explicitly discarded; escape is separately accumulated. The outer event commit only subtracts one from first-passage action.
8. Radius uses exponentially weighted active accumulated slip, not total/wake state or net signed slip.
9. Multiplicity uses radius and a mesh-independent along-front width computed from active unsigned near-tip density; it does not use total or wake content. The v10.2.22 overlay explicitly prevents the ahead-tip MPZ `dx` from acting as the along-front width floor.
10. The sharp-wake PF structural solve does not see analytical tip radius. Radius acts only in local opening/cleavage/emission conversion and is deliberately disabled as a shielding-kernel interpolation coordinate.

At re-initiation the candidate span of the active near-tip unsigned population is only 0.058%; backstress spans 0.056%, radius 0.090%, and signed shielding 0.088%. Thus the orders-of-magnitude total-state diversity is largely orthogonal to the quantities consumed by the current fracture source.
