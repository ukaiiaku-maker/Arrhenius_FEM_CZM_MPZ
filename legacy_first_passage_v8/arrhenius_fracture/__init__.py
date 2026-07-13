"""
Arrhenius phase-field fracture simulation package.

Couples Arrhenius-activated plasticity with AT2 phase-field fracture
for modeling ductile-to-brittle transition in BCC metals (tungsten).
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
from .diagnostics import (
    StepDiagnostics, SimulationHistory,
    save_history, plot_diagnostics, plot_toughness_vs_T,
)

__version__ = '0.1.0'
