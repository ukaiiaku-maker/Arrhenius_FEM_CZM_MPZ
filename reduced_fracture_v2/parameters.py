from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json


@dataclass(frozen=True)
class CanonicalParameters:
    candidate_id: str
    values: dict[str, float]

    @classmethod
    def from_mapping(cls, row: dict[str, str]) -> "CanonicalParameters":
        return cls(str(row["candidate_id"]), {k: float(v) for k, v in row.items() if k != "candidate_id" and _floatable(v)})

    @property
    def fingerprint(self) -> str:
        payload={k:format(v,".17g") for k,v in sorted(self.values.items())}
        return hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(",",":")).encode()).hexdigest()


def _floatable(value: object) -> bool:
    try: float(value); return True
    except (TypeError, ValueError): return False
