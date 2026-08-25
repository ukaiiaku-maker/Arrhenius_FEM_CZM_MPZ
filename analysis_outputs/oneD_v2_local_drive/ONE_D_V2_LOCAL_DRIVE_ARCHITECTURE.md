# One-dimensional V2 local-drive architecture

Status: **ONE_D_FRONT_MODEL_WITH_LOCAL_TENSOR_DRIVE**.

The scalar crack coordinate remains the only geometric evolution coordinate. Global reaction/J/KJ/G values are retained unchanged, but they are explicitly separated from local kinetic fields. `NativeDriveBundle` carries the crack-local basis, a bounded symmetric tensor profile, source-resolved signed emission channels, geometry, quality, and provenance. Scalar KJ alone fails closed for PF signed emission. No singular-field extrapolation is permitted.
