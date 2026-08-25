"""Lossless production-native state factories and dual-lane source shadows.

Each factory constructs one source object and clones it once at initialization.
The source and V2 lanes then evolve independently; state is never copied
between lanes after initialization.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, fields, is_dataclass
from pathlib import Path
from typing import Any, Callable, Mapping
import copy
import hashlib
import json

import numpy as np

from .lifecycle import FEMCZMLifecyclePolicy, PFLifecyclePolicy

STATE_SCHEMA = "oneD_v2_lossless_production_native_state_v1"


def _qualified_name(value: Any) -> str:
    cls = value if isinstance(value, type) else type(value)
    return f"{cls.__module__}.{cls.__qualname__}"


def canonical_source_value(value: Any, seen: set[int] | None = None) -> Any:
    """Return a deterministic, JSON-safe, value-preserving representation."""
    if seen is None:
        seen = set()
    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, (float, np.floating)):
        x = float(value)
        if np.isnan(x): return {"__float__": "nan"}
        if np.isposinf(x): return {"__float__": "+inf"}
        if np.isneginf(x): return {"__float__": "-inf"}
        return x
    if isinstance(value, np.integer): return int(value)
    if isinstance(value, np.ndarray):
        return {"__ndarray__": True, "dtype": str(value.dtype),
                "shape": list(value.shape), "values": value.tolist()}
    if isinstance(value, np.random.Generator):
        return {"__generator__": _qualified_name(value.bit_generator),
                "state": canonical_source_value(value.bit_generator.state, seen)}
    if isinstance(value, Path): return {"__path__": str(value)}
    if isinstance(value, Mapping):
        return {str(k): canonical_source_value(v, seen)
                for k, v in sorted(value.items(), key=lambda item: str(item[0]))}
    if isinstance(value, (list, tuple)):
        return [canonical_source_value(v, seen) for v in value]
    if isinstance(value, (set, frozenset)):
        items = [canonical_source_value(v, seen) for v in value]
        return sorted(items, key=lambda item: json.dumps(item, sort_keys=True, default=str))
    if callable(value):
        return {"__callable__": getattr(value, "__qualname__", repr(value))}
    identity = id(value)
    if identity in seen: return {"__cycle__": _qualified_name(value)}
    seen.add(identity)
    try:
        if is_dataclass(value):
            payload = {f.name: canonical_source_value(getattr(value, f.name), seen)
                       for f in fields(value)}
        elif hasattr(value, "state_dict") and callable(value.state_dict):
            payload = canonical_source_value(value.state_dict(), seen)
        elif hasattr(value, "__dict__"):
            payload = {str(k): canonical_source_value(v, seen)
                       for k, v in sorted(vars(value).items())}
        else: payload = {"repr": repr(value)}
        return {"__source_type__": _qualified_name(value), "fields": payload}
    finally:
        seen.remove(identity)


def scientific_fingerprint(value: Any) -> str:
    encoded = json.dumps(canonical_source_value(value), sort_keys=True,
                         separators=(",", ":"), allow_nan=False).encode()
    return hashlib.sha256(encoded).hexdigest()


def compact_runtime_state(engine: Any) -> Mapping[str, Any]:
    """Bounded exact-comparison view for long predictive trajectories.

    Append-only diagnostic histories are represented by length and terminal
    value; immutable source configuration was already checked by the lossless
    factory export.  All evolving scalars, arrays, process-zone fields, RNG
    states, thresholds, and transaction counters remain value-preserving.
    """
    def compact(value: Any, key: str = "", seen: set[int] | None = None) -> Any:
        if seen is None: seen=set()
        low=key.lower()
        if isinstance(value,(list,tuple)) and ("history" in low or "records" in low):
            return {"length":len(value),"last":None if not value else compact(value[-1],key+".last",seen)}
        if value is None or isinstance(value,(str,bool,int,float,np.integer,np.floating,np.ndarray,np.random.Generator,Path)):
            return canonical_source_value(value)
        if isinstance(value,Mapping):
            return {str(k):compact(v,str(k),seen) for k,v in sorted(value.items(),key=lambda x:str(x[0]))}
        if isinstance(value,(list,tuple)): return [compact(v,key,seen) for v in value]
        ident=id(value)
        if ident in seen: return {"__cycle__":_qualified_name(value)}
        if hasattr(value,"state_dict") and callable(value.state_dict):
            return compact(value.state_dict(),key,seen)
        if hasattr(value,"__dict__"):
            seen.add(ident)
            try:
                return {"__source_type__":_qualified_name(value),"fields":{
                    str(k):compact(v,str(k),seen) for k,v in sorted(vars(value).items())
                    if not callable(v)}}
            finally: seen.remove(ident)
        return canonical_source_value(value)
    excluded={"cb","eb","f","manifest","mpz_config"}
    return {str(k):compact(v,str(k)) for k,v in sorted(vars(engine).items())
            if k not in excluded and not callable(v)}


def _source_capture(engine: Any) -> Mapping[str, Any]:
    capture = getattr(engine, "_capture_state", None)
    transaction = capture() if callable(capture) else copy.deepcopy(vars(engine))
    return {"engine_type": _qualified_name(engine), "transaction_state": transaction,
            "barrier_cleavage": getattr(engine, "cb", None),
            "barrier_emission": getattr(engine, "eb", None),
            "front_config": getattr(engine, "f", None),
            "process_zone_config": getattr(engine, "mpz_config", None),
            "process_zone_state": getattr(engine, "mpz", getattr(engine, "mpz_state", None))}


@dataclass(frozen=True)
class TypedNativeState:
    backend: str
    schema: str
    source_engine_type: str
    source_owned_fields: Mapping[str, Any]
    fingerprint: str

    @classmethod
    def export(cls, backend: str, engine: Any) -> "TypedNativeState":
        payload = canonical_source_value(_source_capture(engine))
        return cls(str(backend), STATE_SCHEMA, _qualified_name(engine), payload,
                   scientific_fingerprint(payload))


@dataclass(frozen=True)
class NativeStateFactoryProduct:
    backend: str
    production_object: Any
    v2_object: Any
    production_state: TypedNativeState
    v2_state: TypedNativeState
    initial_state_fingerprint: str
    independent_lane_objects: bool
    future_state_imports: int = 0


class _NativeStateFactory:
    backend = "UNSPECIFIED"
    def __init__(self, source_constructor: Callable[[], Any]): self.source_constructor = source_constructor
    def instantiate(self) -> NativeStateFactoryProduct:
        source = self.source_constructor(); reduced = copy.deepcopy(source)
        a = TypedNativeState.export(self.backend, source); b = TypedNativeState.export(self.backend, reduced)
        if a.source_owned_fields != b.source_owned_fields:
            raise RuntimeError(f"{self.backend} initialization clone is not lossless")
        if source is reduced: raise RuntimeError("dual lanes must own independent objects")
        return NativeStateFactoryProduct(self.backend, source, reduced, a, b,
                                         a.fingerprint, True, 0)


class PFNativeStateFactory(_NativeStateFactory): backend = "PF"
class FEMCZMNativeStateFactory(_NativeStateFactory): backend = "FEMCZM"


@dataclass(frozen=True)
class DualLaneStep:
    source_result: Mapping[str, Any]
    v2_result: Mapping[str, Any]
    source_state: TypedNativeState
    v2_state: TypedNativeState
    exact_state_match: bool
    maximum_numeric_absolute_error: float
    field_status: str


def _maximum_numeric_error(left: Any, right: Any) -> float:
    errors: list[float] = []
    def visit(a: Any, b: Any) -> None:
        if isinstance(a, Mapping) and isinstance(b, Mapping):
            for key in set(a) & set(b): visit(a[key], b[key])
        elif isinstance(a, list) and isinstance(b, list):
            for x, y in zip(a, b): visit(x, y)
        elif isinstance(a, (int, float)) and isinstance(b, (int, float)):
            errors.append(abs(float(a)-float(b)))
    visit(left, right); return max(errors, default=0.0)


class ExactSourceDualLane:
    """Direct source lane versus V2 lifecycle lane, with no resynchronization."""
    def __init__(self, product: NativeStateFactoryProduct):
        self.product = product; self.backend = product.backend
        self.source = product.production_object; self.v2 = product.v2_object
        self.future_state_imports = 0
        self.lifecycle = PFLifecyclePolicy() if self.backend == "PF" else FEMCZMLifecyclePolicy()
    def advance(self, K_cleave: float, K_emit: float, temperature_K: float,
                dt_s: float, *, full_state: bool = True) -> DualLaneStep:
        if self.backend == "PF":
            source_result = self.source.step(K_cleave, temperature_K, dt_s)
            v2_result = asdict(self.lifecycle.step(self.v2, K_cleave, temperature_K, dt_s))
        else:
            source_result = self.source.step_drives(K_cleave, K_emit, temperature_K, dt_s)
            v2_result = asdict(self.lifecycle.step_engine(self.v2, K_cleave, temperature_K, dt_s))
        if full_state:
            a = TypedNativeState.export(self.backend, self.source); b = TypedNativeState.export(self.backend, self.v2)
        else:
            av=compact_runtime_state(self.source);bv=compact_runtime_state(self.v2)
            a=TypedNativeState(self.backend,"oneD_v2_compact_runtime_state_v1",_qualified_name(self.source),av,scientific_fingerprint(av))
            b=TypedNativeState(self.backend,"oneD_v2_compact_runtime_state_v1",_qualified_name(self.v2),bv,scientific_fingerprint(bv))
        exact = a.source_owned_fields == b.source_owned_fields
        return DualLaneStep(canonical_source_value(source_result), canonical_source_value(v2_result), a, b,
                            exact, _maximum_numeric_error(a.source_owned_fields, b.source_owned_fields),
                            "EXACT_SOURCE_MATCH" if exact else "IMPLEMENTATION_MISMATCH")


def compare_source_fields(source_state: TypedNativeState, v2_state: TypedNativeState) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    def visit(path: str, a: Any, b: Any) -> None:
        if isinstance(a, Mapping) and isinstance(b, Mapping):
            for key in sorted(set(a) | set(b)):
                child = f"{path}.{key}" if path else str(key)
                if key not in a or key not in b: rows.append({"field": child, "status": "IMPLEMENTATION_MISMATCH"})
                else: visit(child, a[key], b[key])
        elif isinstance(a, list) and isinstance(b, list):
            if len(a) != len(b): rows.append({"field": path, "status": "IMPLEMENTATION_MISMATCH"})
            else:
                for i, (x, y) in enumerate(zip(a, b)): visit(f"{path}[{i}]", x, y)
        else:
            rows.append({"field": path, "status": "EXACT_SOURCE_MATCH" if a == b else "IMPLEMENTATION_MISMATCH",
                         "source_value": a, "v2_value": b})
    visit("", source_state.source_owned_fields, v2_state.source_owned_fields); return rows
