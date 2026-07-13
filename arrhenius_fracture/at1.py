"""
AT1 phase-field fracture model
==============================

A dedicated module for the AT1 variant of the variational (Francfort-Marigo /
Bourdin) phase-field fracture model, separated from the AT2 implementation in
``phase_field.py`` so the two degradation/dissipation choices can be developed
and compared cleanly.

Why AT1 (vs AT2)
----------------
AT2 uses a quadratic local dissipation ``w(d) = d^2``; its damage grows from
*zero* applied load (no elastic threshold), which both is unphysical and is a
direct contributor to diffuse damage smearing.  AT1 uses a *linear* local
dissipation ``w(d) = d``.  The nonzero derivative of ``w`` at ``d=0`` produces a
genuine elastic threshold: no damage below a critical stress, a real elastic
stage, and sharper cracks.

Functional and equations
-------------------------
    E[u, d] = integral  (1 - d)^2 * psi_e(u)
                       + (3 Gc / 8) * ( d / ell + ell * |grad d|^2 )   dV

Stationarity in ``d`` (with history field ``H = max_t psi_e`` for irreversible
driving) gives the linear complementarity / variational-inequality problem

    2 diag(H) Md d  +  (3 Gc ell / 4) Kd d   =   2 Md H  -  (3 Gc / (8 ell)) Md 1
    subject to    d_prev <= d <= 1         (irreversibility + saturation)

solved here by a primal active-set iteration on the lower bound.  The constant
force term ``-(3 Gc/(8 ell)) Md 1`` (from the linear local term) is what makes
``d = d_prev`` the exact solution below threshold -- i.e. the elastic stage.

Closed-form properties (used as unit checks and as diagnostics)
---------------------------------------------------------------
Critical elastic energy density for damage onset:

    H_c = 3 Gc / (16 ell)

Critical (peak) stress -- the material STRENGTH implied by (Gc, ell):

    sigma_c = sqrt( 3 E Gc / (8 ell) )

Homogeneous 1-D damage above threshold (uniform field, |grad d| = 0):

    d(H) = 1 - H_c / H        for H >= H_c,    d = 0 otherwise

The implementation reproduces ``d(H) = 1 - H_c/H`` exactly (verified).

-------------------------------------------------------------------------------
NOTES / KNOWN LIMITATION -- the strength-ell coupling  (and PF-CZM as a later fix)
-------------------------------------------------------------------------------
AT1 ties the material strength to the regularization length:
``sigma_c = sqrt(3 E Gc / (8 ell))``.  On a mesh that resolves ``ell`` (need
``ell >~ 2 h``), the achievable ``ell`` sets ``sigma_c``, and for the W SENT case
this lands ``sigma_c`` right at the far-field stress that drives the crack:

    ell ~ 3 um  (finest on the adaptive ~1.5 um tip mesh)
        -> sigma_c ~ sqrt(3 * 4.5e11 * 1.39 / (8 * 3e-6)) ~ 280 MPa
    K = 13 MPa.sqrt(m) on this SENT  -> sigma_far ~ 290 MPa

so at the crack-driving load the WHOLE specimen reaches its own strength and
damages in bulk instead of localizing a crack.  Clean localization needs
``sigma_c`` to clear ``sigma_far`` by ~2-3x, i.e. ``ell <~ 0.7 um``, i.e. a
sub-micron mesh ALONG THE WHOLE CRACK PATH (tip-only refinement is not enough --
the crack goes blocky once it leaves the refined zone).

==> Possible later upgrade: **PF-CZM (Wu, 2017)** phase-field cohesive-zone model.
    It makes the strength ``sigma_c`` (a.k.a. f_t) an INDEPENDENT material input,
    decoupled from ``ell`` (which becomes a pure numerical regularizer).  That is
    the principled fix for the coupling above: set sigma_c to the physical W value
    (GPa) and keep ``ell`` only for resolution, so damage localizes at the tip on
    an affordable mesh.  Switching is a change to the degradation function g(d)
    and the dissipation w(d)/c_w (the active-set solver scaffolding here is
    largely reusable).  We start with AT1; revisit PF-CZM if the strength-ell
    coupling proves limiting.

Coupling to the Arrhenius kinetics
----------------------------------
Standard phase-field is rate-independent (quasi-static energy minimization).
The Arrhenius / first-passage MOBILITY coupling (the rate-/T-dependence that
spans DBTT-to-ceramic) enters through ``drive_multiplier`` (a tip-localized
amplification of the elastic driving energy ``H``).  IMPORTANT: that multiplier
must be tip-localized -- a global ``M_fract > 1`` over the bulk pushes the whole
field over the AT1 threshold and nucleates damage everywhere.  The default here
is pure-variational (multiplier off); the kinetic coupling is opt-in.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from .mesh import TriMesh


# --------------------------------------------------------------------------- #
#  Closed-form AT1 diagnostics (also used as unit checks)
# --------------------------------------------------------------------------- #
def at1_critical_energy_density(Gc, ell: float):
    """H_c = 3 Gc / (16 ell): elastic-energy-density threshold for damage onset."""
    return 3.0 * np.asarray(Gc, dtype=float) / (16.0 * ell)


def at1_critical_stress(E: float, Gc, ell: float):
    """sigma_c = sqrt(3 E Gc / (8 ell)): the material strength implied by (E, Gc, ell)."""
    return np.sqrt(3.0 * E * np.asarray(Gc, dtype=float) / (8.0 * ell))


def at1_surface_energy(mesh: TriMesh, d: np.ndarray, ell: float, Gc) -> float:
    """AT1 crack surface energy   E = integral Gc * (3/8)(d/ell + ell|grad d|^2) dV.

    Gc may be scalar or per-node.  Mirrors ``phase_field.at2_surface_energy`` but
    with the LINEAR local term (integral of d, not d^2)."""
    Gc_arr = np.broadcast_to(np.asarray(Gc, dtype=float), (mesh.nn,))
    E = 0.0
    for e in range(mesh.ne):
        conn = mesh.elems[e]
        de = d[conn]
        A = mesh.area_e[e]
        Gc_e = float(np.mean(Gc_arr[conn]))
        # integral of a linear field over a triangle = A/3 * sum(nodal values)
        int_d = A / 3.0 * (de[0] + de[1] + de[2])
        grad_d = mesh.dNdx_e[e] @ de
        int_g2 = A * (grad_d[0] ** 2 + grad_d[1] ** 2)
        E += Gc_e * (3.0 / 8.0) * ((1.0 / ell) * int_d + ell * int_g2)
    return E


# --------------------------------------------------------------------------- #
#  AT1 damage update (bound-constrained variational solve)
# --------------------------------------------------------------------------- #
def update_phase_field_at1(
    d, Hhist, psi_e_node, Md, Kd, notch_nodes, Gc_eff, ell, Gamma0, dt,
    drive_multiplier=None, cohesive_gate=None, use_kinetic_drive=False,
    max_damage_increment=None, max_active_set_iter=15,
):
    """Clean AT1 phase-field update (linear local dissipation w(d)=d).

    Solves   2 diag(H) Md d + (3 Gc ell/4) Kd d = 2 Md H - (3 Gc/(8 ell)) Md 1
    with bound constraints  d_prev <= d <= 1  via a primal active-set iteration
    on the lower bound.  Below ``H_c = 3 Gc/(16 ell)`` the unconstrained solution
    is negative, so the active set pins ``d = d_prev`` (= 0 in virgin material):
    the elastic threshold, as a variational consequence of w(d)=d -- NOT a gate.

    No source-mode caps/hazards/floors.  ``drive_multiplier`` (tip-localized, see
    module notes) is the optional (A) kinetic coupling; default pure-variational.

    Returns (d_new, Hhist).
    """
    nn = len(d)
    Hhist = np.maximum(Hhist, psi_e_node)
    Gc = np.maximum(np.broadcast_to(np.asarray(Gc_eff, float), (nn,)).copy(), 1e-30)

    Heff = Hhist
    if drive_multiplier is not None:
        dm = np.clip(np.broadcast_to(np.asarray(drive_multiplier, float), (nn,)), 0.0, 1e3)
        Heff = dm * Hhist
    if cohesive_gate is not None:
        Heff = Heff * np.clip(np.broadcast_to(np.asarray(cohesive_gate, float), (nn,)), 0.0, 1.0)

    A = (sparse.diags(3.0 * Gc * ell / 4.0) @ Kd) + 2.0 * (sparse.diags(Heff) @ Md)
    b = 2.0 * (Md @ Heff) - (Md @ (3.0 * Gc / (8.0 * ell)))

    d_lb = d.copy(); d_lb[notch_nodes] = 1.0           # irreversibility lower bound
    d_out = np.maximum(d.copy(), 0.0); d_out[notch_nodes] = 1.0
    active_lo = np.zeros(nn, dtype=bool)
    active_lo[notch_nodes] = True                      # notch pinned at lb (= 1)

    for _ in range(max_active_set_iter):
        solve = ~active_lo
        if not solve.any():
            break
        d_fix = np.where(active_lo, d_lb, 0.0)
        rhs = b - A @ d_fix
        d_new = d_out.copy()
        d_new[solve] = spsolve(A[np.ix_(solve, solve)], rhs[solve])
        d_new[active_lo] = d_lb[active_lo]
        viol = solve & (d_new < d_lb - 1e-12)          # nodes wanting to drop below lb
        if not viol.any():
            d_out = d_new
            break
        active_lo |= viol
        d_out = d_new

    d_out = np.clip(d_out, d_lb, 1.0)

    if (max_damage_increment is not None and np.isfinite(max_damage_increment)
            and max_damage_increment > 0):
        d_out = np.minimum(d_out, d + max_damage_increment)
        d_out = np.clip(d_out, d_lb, 1.0)

    return d_out, Hhist
