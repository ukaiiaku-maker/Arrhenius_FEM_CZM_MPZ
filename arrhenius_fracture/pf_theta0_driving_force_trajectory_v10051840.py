"""Audited PF v10.4.1 driving-force trajectory reader for theta=0 parity.

Reads the reference `steps_1000K.csv` step table and exposes an immutable,
provenance-stamped `J_target(t)` / `KJ_target(t)` history that the
transactional loading controller (Phase 2) consumes. This module performs no
FEM state mutation and no hazard/RNG interaction; it is a pure reader.

Physical time is reconstructed from the recorded per-row `dt_cur_s` interval,
never from solver step count alone, so that a target time genuinely
corresponds to a physical time offset in the PF reference run.

The "prefracture" slice is the strict prefix of rows before the PF reference
first commits nonzero `crack_extension_m`. Per FEM_CZM_HANDOFF.md section 6,
a stochastic crossing event must never be assigned to an interval endpoint,
so the prefracture slice deliberately excludes the first-passage row itself;
callers that need the first-passage target values read them from the
`first_passage_*` fields instead.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd

MODEL_ID = "pf_theta0_driving_force_trajectory_v10_0_5_18_4_0"

REQUIRED_COLUMNS = (
    "step",
    "dt_cur_s",
    "J_effective_direct_J_per_m2",
    "J_signed_direct_J_per_m2",
    "KJ_Pa_sqrtm",
    "B",
    "N_em",
    "crack_extension_m",
)


@dataclass(frozen=True, eq=False)
class PFDrivingForceTrajectory:
    schema: str
    source_path: str
    source_sha256: str
    n_rows: int
    step: np.ndarray
    physical_time_s: np.ndarray
    dt_s: np.ndarray
    J_target_J_per_m2: np.ndarray
    J_signed_J_per_m2: np.ndarray
    KJ_target_Pa_sqrtm: np.ndarray
    B: np.ndarray
    N_em: np.ndarray
    crack_extension_m: np.ndarray
    W_bulk_plastic_cumulative_J_per_m: np.ndarray
    first_passage_row_index: int | None
    first_passage_step: float | None
    first_passage_time_s: float | None
    first_passage_J_J_per_m2: float | None
    first_passage_KJ_Pa_sqrtm: float | None


def _sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _column(frame: pd.DataFrame, name: str) -> np.ndarray:
    return pd.to_numeric(frame[name], errors="raise").to_numpy(dtype=float)


def load_pf_theta0_driving_force_trajectory(
    steps_csv_path: Path | str,
) -> PFDrivingForceTrajectory:
    path = Path(steps_csv_path).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"PF steps table not found: {path}")

    frame = pd.read_csv(path)
    if frame.empty:
        raise ValueError(f"empty PF steps table: {path}")

    missing = [name for name in REQUIRED_COLUMNS if name not in frame.columns]
    if missing:
        raise ValueError(f"{path} lacks required columns: {sorted(missing)}")

    step = _column(frame, "step")
    dt_s = _column(frame, "dt_cur_s")
    if np.any(dt_s <= 0.0):
        raise ValueError(f"{path} contains a non-positive dt_cur_s interval")
    if np.any(np.diff(step) <= 0.0):
        raise ValueError(f"{path} step column is not strictly increasing")

    physical_time_s = np.cumsum(dt_s)

    J_target = _column(frame, "J_effective_direct_J_per_m2")
    J_signed = _column(frame, "J_signed_direct_J_per_m2")
    KJ_target = _column(frame, "KJ_Pa_sqrtm")
    B = _column(frame, "B")
    N_em = _column(frame, "N_em")
    crack_extension_m = _column(frame, "crack_extension_m")
    if "W_bulk_plastic_cumulative_J_per_m" in frame.columns:
        W_bulk_plastic = _column(frame, "W_bulk_plastic_cumulative_J_per_m")
    else:
        W_bulk_plastic = np.zeros_like(step)

    first_passage_row_index: int | None = None
    nonzero_extension = np.flatnonzero(crack_extension_m > 0.0)
    if nonzero_extension.size:
        first_passage_row_index = int(nonzero_extension[0])

    if first_passage_row_index is None:
        first_passage_step = None
        first_passage_time_s = None
        first_passage_J = None
        first_passage_KJ = None
    else:
        i = first_passage_row_index
        first_passage_step = float(step[i])
        first_passage_time_s = float(physical_time_s[i])
        first_passage_J = float(J_target[i])
        first_passage_KJ = float(KJ_target[i])

    return PFDrivingForceTrajectory(
        schema=MODEL_ID,
        source_path=str(path),
        source_sha256=_sha256_of(path),
        n_rows=int(len(frame)),
        step=step,
        physical_time_s=physical_time_s,
        dt_s=dt_s,
        J_target_J_per_m2=J_target,
        J_signed_J_per_m2=J_signed,
        KJ_target_Pa_sqrtm=KJ_target,
        B=B,
        N_em=N_em,
        crack_extension_m=crack_extension_m,
        W_bulk_plastic_cumulative_J_per_m=W_bulk_plastic,
        first_passage_row_index=first_passage_row_index,
        first_passage_step=first_passage_step,
        first_passage_time_s=first_passage_time_s,
        first_passage_J_J_per_m2=first_passage_J,
        first_passage_KJ_Pa_sqrtm=first_passage_KJ,
    )


def prefracture_slice(trajectory: PFDrivingForceTrajectory) -> PFDrivingForceTrajectory:
    """Return the strict prefix of rows before the PF reference's first
    nonzero crack_extension_m. The first-passage row itself is excluded: the
    controller must localize a crossing event inside an interval, never
    treat the interval endpoint as the event (FEM_CZM_HANDOFF.md section 6).
    """
    end = (
        trajectory.n_rows
        if trajectory.first_passage_row_index is None
        else trajectory.first_passage_row_index
    )
    return PFDrivingForceTrajectory(
        schema=trajectory.schema,
        source_path=trajectory.source_path,
        source_sha256=trajectory.source_sha256,
        n_rows=end,
        step=trajectory.step[:end],
        physical_time_s=trajectory.physical_time_s[:end],
        dt_s=trajectory.dt_s[:end],
        J_target_J_per_m2=trajectory.J_target_J_per_m2[:end],
        J_signed_J_per_m2=trajectory.J_signed_J_per_m2[:end],
        KJ_target_Pa_sqrtm=trajectory.KJ_target_Pa_sqrtm[:end],
        B=trajectory.B[:end],
        N_em=trajectory.N_em[:end],
        crack_extension_m=trajectory.crack_extension_m[:end],
        W_bulk_plastic_cumulative_J_per_m=trajectory.W_bulk_plastic_cumulative_J_per_m[:end],
        first_passage_row_index=trajectory.first_passage_row_index,
        first_passage_step=trajectory.first_passage_step,
        first_passage_time_s=trajectory.first_passage_time_s,
        first_passage_J_J_per_m2=trajectory.first_passage_J_J_per_m2,
        first_passage_KJ_Pa_sqrtm=trajectory.first_passage_KJ_Pa_sqrtm,
    )


_CHANNELS = {
    "J": "J_target_J_per_m2",
    "KJ": "KJ_target_Pa_sqrtm",
}


def interpolate_target(
    trajectory: PFDrivingForceTrajectory,
    physical_time_s: float,
    *,
    channel: str = "KJ",
) -> float:
    """Bounded linear interpolation of J or KJ at a physical time.

    Raises ValueError outside [time[0], time[-1]]: extrapolation of the PF
    reference target is not permitted, matching the no-extrapolation
    contract already used for the PF signed shielding kernel.
    """
    if channel not in _CHANNELS:
        raise ValueError(f"unknown channel {channel!r}; expected one of {sorted(_CHANNELS)}")
    times = trajectory.physical_time_s
    lo, hi = float(times[0]), float(times[-1])
    if physical_time_s < lo or physical_time_s > hi:
        raise ValueError(
            f"requested physical_time_s={physical_time_s} outside PF reference "
            f"target range [{lo}, {hi}]; extrapolation is not permitted"
        )
    values = getattr(trajectory, _CHANNELS[channel])
    return float(np.interp(physical_time_s, times, values))


__all__ = [
    "MODEL_ID",
    "REQUIRED_COLUMNS",
    "PFDrivingForceTrajectory",
    "load_pf_theta0_driving_force_trajectory",
    "prefracture_slice",
    "interpolate_target",
]
