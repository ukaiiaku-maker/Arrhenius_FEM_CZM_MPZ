# V2 exact mechanics-oracle architecture

Status: **IMPLEMENTED, FAIL-CLOSED**.

`BackendElasticFieldOracle` returns a candidate-independent `NormalizedElasticFieldSnapshot`; `BackendSourceProbeOperator` applies the current physical process-zone state and returns `NativeDriveBundle`. Field keys include source/factory, full geometry, mesh/topology, crack basis and normalization. Probe keys separately include snapshot, radius, width, area, source geometry/positions, slip-system identity and convention. Backstress, shielding, populations and candidate identity are forbidden from field snapshots. `OptionalRegimeAwareSurrogate` is locked while Level 3 is not qualified.
