# PF tensor interpolation failure audit

Status: **SPARSE EXTENSION-ONLY INTERPOLATION REJECTED**.

All 23 sign failures were re-evaluated against the realized wake, front basis,
probe support, probe weights, reliability and categorical signed-system state.
They are not one undifferentiated numerical error: `{"FRONT_DIRECTION_TRANSITION": 9, "NEAR_ZERO_COMPONENT_SIGN_ARTIFACT": 5, "PROBE_ELEMENT_SELECTION_TRANSITION": 9}`.
The CSV retains all simultaneous cause flags and the left/target/right regime
keys. A reliability, signed-system, or sign transition is a hard failure. A
near-zero component is reported as such and is never promoted to a reliable
continuous sign crossing.
