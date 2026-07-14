"""Arrhenius hazard-based fracture and fatigue simulation package.

Version 0.9 adds a conservative moving one-dimensional crack-tip process zone
while retaining the anisotropic, multifront, branching, coalescence, cyclic
mechanics, mixed-mode, sharp-wake, and adaptive-CZM production architecture.
The frozen scalar v8 closure remains selectable as ``legacy_scalar``.
"""

from .config import (
    SimulationConfig, GeometryConfig, MeshConfig, ElasticProperties,
    PlasticityBarrier, FractureBarrier, LoadingConfig, PhaseFieldConfig,
    DislocationConfig, TipMemoryConfig, CohesiveConfig, HazardConfig,
    JIntegralConfig, DiagnosticsConfig,
    make_dbtt_config, make_ceramic_config, make_cohesive_dbtt_config,
    make_emergent_config,
    KB, EV_TO_J,
)
from .mesh import TriMesh, BoundaryData, make_tri_mesh, make_boundary_data
from .materials import PlasticityModel, FractureModel
from .fem import (
    plane_strain_D, assemble_mechanics, solve_dirichlet,
    project_gp_to_nodes, assemble_pf_matrices,
)
from .plasticity import update_plasticity
from .phase_field import (
    update_phase_field, at2_surface_energy, crack_front_mask,
    update_tip_memory, compute_fracture_amplification, cohesive_gate,
)
from .j_integral import compute_J_integral, find_crack_tip, compute_crack_advance

from .moving_process_zone import MovingProcessZoneConfig, MovingProcessZoneState
from .diagnostics import (
    StepDiagnostics, SimulationHistory,
    save_history, plot_diagnostics, plot_toughness_vs_T,
)

# Preserve the existing import surface while activating the v9.4 signed
# detailed-balance completion everywhere the v9.3 module is imported lazily.
from . import emission_derived_plasticity as _pt_v93
from .emission_derived_plasticity_v94 import (
    EmissionDerivedPeierlsTaylorModel as _PTModelV94,
)


def _exact_pt_series_rate(rate_a, rate_b):
    """Harmonic sequential rate with exact absorbing zero bottlenecks."""
    import numpy as np

    a, b = np.broadcast_arrays(
        np.maximum(np.asarray(rate_a, dtype=float), 0.0),
        np.maximum(np.asarray(rate_b, dtype=float), 0.0),
    )
    out = np.zeros_like(a, dtype=float)
    active = (a > 0.0) & (b > 0.0)
    np.divide(a * b, a + b, out=out, where=active)
    return out


_PTModelV94.series_rate = staticmethod(_exact_pt_series_rate)
_pt_v93.EmissionDerivedPeierlsTaylorModel = _PTModelV94

# v9.5 keeps the v9.4 constitutive closure but evaluates forest density and
# effective stress locally in every moving-process-zone bin.  Patch the module
# attribute before mpz_front_engine imports it, while preserving public imports.
from . import moving_process_zone as _mpz_v94
from .moving_process_zone_v95 import MovingProcessZoneState as _MPZStateV95
_mpz_v94.MovingProcessZoneState = _MPZStateV95
MovingProcessZoneState = _MPZStateV95

__version__ = '0.9.5'
