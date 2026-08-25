# PF tensor-map qualification

Status: **SOURCE NODES QUALIFIED; SPARSE INTERPOLATION UNQUALIFIED**.

The archives preserve derived signed shears and factors, but not raw tensors or probe weights. Two bounded default-off observer diagnostics therefore captured the missing fields. Peak seed 8666 stopped after five accepted geometry events at 14.508 µm; DBTT seed 1008666 stopped after three at 10.475 µm. Across the physical CSV/JSON artifacts, OFF and ON hashes are byte-identical. The observer recorded 216 Peak and 247 DBTT tensor evaluations.

A separate deterministic candidate-independent replay used the production plane-strain graded mesh, finite-width binary stiffness-killed wake, theta=0 cubic elasticity, the corrected front-direction selector, and exact source probe. It stores 376 raw selected-element profile rows at the pinned 0–1000 µm map nodes. Every source probe was reliable; reaction, all tensor components, and resolved signed shear scale exactly with opening (maximum error 0).

The sparse node set is nevertheless inadequate for forward interpolation. Leave-one-node-out validation produced 23 component/projection sign mismatches and a maximum relative error of 2.2033e4. This is the discrete wake/mesh topology signal, not a smooth continuum profile. Therefore the parquet and node summary are qualified source-node evidence, but `PFTensorDriveProvider` rejects the map for general forward use. No compression or extrapolation is used.
