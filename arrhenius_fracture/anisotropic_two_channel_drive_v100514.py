"""PF v10.2.22-style signed two-channel drives from FEM tensor probes.

The opening tensor is sampled ahead of the crack front.  Each reduced BCC slip
trace is sampled along its own forward ray, matching the reference PF geometry
rather than resolving both channels from one common tensor.
"""
from __future__ import annotations

import copy
from typing import Any, Callable, Sequence

import numpy as np

from .crystal import bcc_slip_traces
from . import mixed_mode_first_passage_v8 as _mm

MODEL_ID = "FEM_CZM_two_reduced_signed_slip_channels_v10_0_5_14"


def _unit(vector: Sequence[float]) -> np.ndarray:
    value = np.asarray(vector, dtype=float).reshape(2)
    norm = float(np.linalg.norm(value))
    if not np.isfinite(norm) or norm <= 1.0e-30:
        raise ValueError("drive direction must be finite and nonzero")
    return value / norm


def resolve_two_channel_drive_from_tensors(
    opening_tensor: np.ndarray,
    channel_tensors: Sequence[np.ndarray],
    crack_direction: Sequence[float],
    crystal_theta_deg: float,
    schmid_reference: float = 0.5,
) -> dict[str, Any]:
    """Resolve signed shear and nonnegative drive factors for two slip traces."""
    opening = np.asarray(opening_tensor, dtype=float).reshape(2, 2)
    if not np.all(np.isfinite(opening)):
        raise ValueError("opening tensor must be finite")
    traces = bcc_slip_traces(float(crystal_theta_deg))
    if len(traces) != 2:
        raise RuntimeError(
            f"v10.0.5.14 requires two reduced BCC traces; got {len(traces)}"
        )
    if len(channel_tensors) != len(traces):
        raise ValueError("one FEM tensor is required for each reduced slip channel")

    forward = _unit(crack_direction)
    normal = np.array([-forward[1], forward[0]], dtype=float)
    sigma_nn = float(normal @ opening @ normal)
    eig = np.linalg.eigvalsh(opening)
    sigma_amplitude = max(float(eig[-1]), max(sigma_nn, 0.0), 1.0)
    reference = max(abs(float(schmid_reference)), 1.0e-12)

    names: list[str] = []
    signed_tau: list[float] = []
    factors: list[float] = []
    trace_directions: list[list[float]] = []
    trace_normals: list[list[float]] = []
    for trace, raw_tensor in zip(traces, channel_tensors):
        tensor = np.asarray(raw_tensor, dtype=float).reshape(2, 2)
        if not np.all(np.isfinite(tensor)):
            raise ValueError("channel tensor must be finite")
        t = _unit(trace["t"])
        n = _unit(trace["n"])
        tau = float(t @ tensor @ n)
        names.append(str(trace["name"]))
        signed_tau.append(tau)
        factors.append(abs(tau) / (reference * sigma_amplitude))
        trace_directions.append(t.tolist())
        trace_normals.append(n.tolist())

    return {
        "two_channel_model": MODEL_ID,
        "two_channel_drive_reliable": True,
        "two_channel_names": names,
        "two_channel_trace_directions": trace_directions,
        "two_channel_trace_normals": trace_normals,
        "two_channel_tau_signed_Pa": signed_tau,
        "two_channel_drive_factors": factors,
        "two_channel_sigma_nn_Pa": sigma_nn,
        "two_channel_sigma_amplitude_Pa": sigma_amplitude,
        "two_channel_schmid_reference": reference,
        "two_channel_factors_normalized": False,
        "two_channel_factors_clipped": False,
    }


def resolve_two_channel_drive(
    stress_tensor: np.ndarray,
    sigma_opening_Pa: float,
    crystal_theta_deg: float,
    schmid_reference: float = 0.5,
) -> dict[str, Any]:
    """Compatibility helper for unit tests using one prescribed tensor.

    Production uses :func:`resolve_two_channel_drive_from_tensors` with separate
    channel-ray probes.  This helper supplies the same tensor for both channels.
    """
    tensor = np.asarray(stress_tensor, dtype=float).reshape(2, 2)
    opening = tensor.copy()
    opening[1, 1] = float(sigma_opening_Pa)
    return resolve_two_channel_drive_from_tensors(
        opening,
        [tensor, tensor],
        [1.0, 0.0],
        crystal_theta_deg,
        schmid_reference,
    )


def _argument(args: tuple[Any, ...], kwargs: dict[str, Any], key: str, index: int):
    if key in kwargs:
        return kwargs[key]
    if len(args) <= index:
        raise TypeError(f"missing required J-wrapper argument {key!r}")
    return args[index]


def augmented_j_wrapper_factory(original_factory: Callable) -> Callable:
    """Preserve the existing J/profile wrapper and add two signed channel probes."""

    def factory(original_compute, context):
        base_wrapped = original_factory(original_compute, context)

        def wrapped(*args, **kwargs):
            result = base_wrapped(*args, **kwargs)
            try:
                mesh = _argument(args, kwargs, "mesh", 0)
                sigma_gp = _argument(args, kwargs, "sigma_gp", 2)
                damage = _argument(args, kwargs, "d", 4)
                crack_tip = _argument(args, kwargs, "crack_tip", 5)
                crack_direction = _unit(
                    _argument(args, kwargs, "crack_direction", 6)
                )
                probe_kwargs = dict(
                    radius_m=float(context.probe_radius_m),
                    annulus_half_width=float(context.annulus_half_width),
                    sector_half_angle_deg=float(context.sector_half_angle_deg),
                    damage_cutoff=float(context.damage_cutoff),
                )
                opening_probe = _mm.process_zone_traction_probe(
                    mesh,
                    sigma_gp,
                    damage,
                    crack_tip,
                    crack_direction,
                    **probe_kwargs,
                )
                if not bool(opening_probe.get("reliable", False)):
                    raise RuntimeError("opening tensor probe is unreliable")

                channel_tensors: list[np.ndarray] = []
                channel_elements: list[int] = []
                channel_expansions: list[float] = []
                for trace in bcc_slip_traces(float(context.crystal_theta_deg)):
                    ray = _unit(trace["t"])
                    if float(ray @ crack_direction) < 0.0:
                        ray = -ray
                    probe = _mm.process_zone_traction_probe(
                        mesh,
                        sigma_gp,
                        damage,
                        crack_tip,
                        ray,
                        **probe_kwargs,
                    )
                    if not bool(probe.get("reliable", False)):
                        raise RuntimeError(
                            f"channel tensor probe is unreliable for {trace['name']}"
                        )
                    channel_tensors.append(
                        np.asarray(probe["stress_tensor"], dtype=float)
                    )
                    channel_elements.append(int(probe.get("n_elements", 0)))
                    channel_expansions.append(float(probe.get("expansion", 0.0)))

                drive = resolve_two_channel_drive_from_tensors(
                    np.asarray(opening_probe["stress_tensor"], dtype=float),
                    channel_tensors,
                    crack_direction,
                    float(context.crystal_theta_deg),
                )
                drive.update(
                    {
                        "two_channel_opening_probe_elements": int(
                            opening_probe.get("n_elements", 0)
                        ),
                        "two_channel_opening_probe_expansion": float(
                            opening_probe.get("expansion", 0.0)
                        ),
                        "two_channel_probe_elements": channel_elements,
                        "two_channel_probe_expansion": channel_expansions,
                    }
                )
            except Exception as exc:
                drive = {
                    "two_channel_model": MODEL_ID,
                    "two_channel_drive_reliable": False,
                    "two_channel_names": [],
                    "two_channel_tau_signed_Pa": [],
                    "two_channel_drive_factors": [],
                    "two_channel_failure": f"{type(exc).__name__}: {exc}",
                }
            context.latest.update(copy.deepcopy(drive))
            if context.records:
                context.records[-1].update(copy.deepcopy(drive))
            try:
                result[2].update(copy.deepcopy(drive))
            except Exception:
                pass
            return result

        return wrapped

    return factory


__all__ = [
    "MODEL_ID",
    "resolve_two_channel_drive",
    "resolve_two_channel_drive_from_tensors",
    "augmented_j_wrapper_factory",
]
