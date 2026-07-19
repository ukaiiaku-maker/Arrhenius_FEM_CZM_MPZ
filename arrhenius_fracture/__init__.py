"""Arrhenius hazard-based fracture and fatigue simulation package.

Version 10.0 adds the selectable ``kinetic_campaign_czm`` front state and its
transactional trial-cohesive infrastructure.  The existing ``legacy_scalar``
and ``moving_pz`` implementations remain available and the package-level
moving-PZ default is intentionally unchanged for regression compatibility.
"""

from .config import (
    SimulationConfig, GeometryConfig, MeshConfig, ElasticProperties,
    PlasticityBarrier, FractureBarrier, LoadingConfig, FractureResistanceConfig,
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
    project_gp_to_nodes,
)
from .plasticity import update_plasticity
from .j_integral import compute_J_integral, find_crack_tip, compute_crack_advance

from .moving_process_zone import MovingProcessZoneConfig, MovingProcessZoneState
from .diagnostics import (
    StepDiagnostics, SimulationHistory,
    save_history, plot_diagnostics, plot_toughness_vs_T,
)

# Preserve the existing import surface while activating the uncapped v9.6
# emission-derived Peierls--Taylor closure everywhere the base module is
# imported lazily by bulk FEM, fatigue, or moving-process-zone paths.
from . import emission_derived_plasticity as _pt_base
from .emission_derived_plasticity_v96 import (
    EmissionDerivedPeierlsTaylorModel as _PTModelV96,
)


def _exact_uncapped_pt_series_rate(rate_a, rate_b):
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


_PTModelV96.series_rate = staticmethod(_exact_uncapped_pt_series_rate)
_pt_base.EmissionDerivedPeierlsTaylorModel = _PTModelV96

# v9.5 remains the package-level moving-PZ default.  The new kinetic state is
# selected explicitly by its v10 entry point; it never silently changes legacy
# runs or old parameterizations.
from . import moving_process_zone as _mpz_base
from .moving_process_zone_v95 import MovingProcessZoneState as _MPZStateV95
_mpz_base.MovingProcessZoneState = _MPZStateV95
MovingProcessZoneState = _MPZStateV95

from .pf_equivalent_material_manifest import (
    PFEquivalentMaterialManifest,
    load_material_manifest as load_pf_equivalent_material_manifest,
)
from .kinetic_campaign_czm import (
    KineticCampaignCZMConfig,
    CampaignKineticMPZState,
    CampaignCalibratedCZMFrontEngine,
    DevelopedStateDiagnosticCZMFrontEngine,
)
from .cohesive_trial_state import (
    KineticCZMTransactionSnapshot,
    KineticTrialAdaptiveCZMBackend,
)
from .kinetic_cohesive_stepper import (
    KineticCohesiveStepperConfig,
    KineticCohesiveStepResult,
    KineticCohesiveStepper,
)

__version__ = '10.0.5.7'
