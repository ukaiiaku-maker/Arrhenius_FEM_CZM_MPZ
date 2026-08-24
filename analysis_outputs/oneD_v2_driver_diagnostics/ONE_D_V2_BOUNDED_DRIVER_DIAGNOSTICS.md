# One D V2 Bounded Driver Diagnostics

The authorized diagnostics were provenance-resolved to PF commit `ab627933` and FEM/CZM commit `931bed6`. Before launch, direct production-source inspection found an irreducible PF rollback failure: late geometry veto invokes `restore_geometry_veto`, which intentionally raises because coupled-step replay is unavailable. A forced-veto diagnostic therefore cannot satisfy full-state restoration. No run was launched or extended, and no trajectory result is claimed.
