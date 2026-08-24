# One-dimensional V2 FEM/CZM transaction lineage audit

Ancestry `30b53ff -> 931bed66913afc970117bce900805ecc9b6225f8`: **True**.

| File | Function / boundary | 30b53ff SHA-256 | 931bed6 SHA-256 | Classification |
|---|---|---|---|---|
| `arrhenius_fracture/bulk_pt_live_run2d_composition_v1043.py` | `build_bulk_run_2d / outer transaction begin` | `3a21b36e2efe37d172b0a626b77fa7ccbbe6ee7fd39133009f60076a71e98bea` | `3a21b36e2efe37d172b0a626b77fa7ccbbe6ee7fd39133009f60076a71e98bea` | **IDENTICAL** |
| `arrhenius_fracture/bulk_pt_live_run2d_composition_v1043.py` | `build_bulk_run_2d / mechanics interval begin` | `feffeb9abd5958cb6bfb217d20de5faa2619988c4383952d57827b7a7fbd325d` | `feffeb9abd5958cb6bfb217d20de5faa2619988c4383952d57827b7a7fbd325d` | **IDENTICAL** |
| `arrhenius_fracture/bulk_pt_live_run2d_composition_v1043.py` | `build_bulk_run_2d / accepted mechanics interval` | `cbfe1612044f56e750bf3415f43ea09e888b952a1cc12fd15feb3a966bf5eb60` | `cbfe1612044f56e750bf3415f43ea09e888b952a1cc12fd15feb3a966bf5eb60` | **IDENTICAL** |
| `arrhenius_fracture/bulk_pt_live_run2d_composition_v1043.py` | `build_bulk_run_2d / pre-remap cohesive snapshot` | `875b6a7aa25a0c00536a3036db4a64b1c3eaa92055b8fbdb9d11aa43742c5088` | `875b6a7aa25a0c00536a3036db4a64b1c3eaa92055b8fbdb9d11aa43742c5088` | **IDENTICAL** |
| `arrhenius_fracture/bulk_pt_live_run2d_composition_v1043.py` | `build_bulk_run_2d / overlap topology remap` | `b89261df5589658d7266e74cab70a02420857a0b4983edb8b1704c2d3e922191` | `b89261df5589658d7266e74cab70a02420857a0b4983edb8b1704c2d3e922191` | **IDENTICAL** |
| `arrhenius_fracture/bulk_pt_live_run2d_composition_v1043.py` | `build_bulk_run_2d / deflect geometry veto rollback` | `6c087f554033a8b6f0d802968be46dfc18ae157a868ec174dc377266cf018b15` | `6c087f554033a8b6f0d802968be46dfc18ae157a868ec174dc377266cf018b15` | **IDENTICAL** |
| `arrhenius_fracture/bulk_pt_live_run2d_composition_v1043.py` | `build_bulk_run_2d / authoritative commit` | `53ad2fb3e96fc7b2188343a839ddb3d7b787f6344ec3c4e1a8407fc7dc5e078c` | `296a25b6cbb21c9f2d84dff9f2451e10efb4613835a6708b4c47652bebbec18d` | **NONSEMANTIC_CHANGE** |
| `arrhenius_fracture/bulk_pt_live_run2d_composition_v1043.py` | `remesh_with_failure_capture / tentative remap/failure capture` | `a54b3c5f0326f686a119e84be46249d40b7b92203ee9907d78fa1eb86426a984` | `74727c01a140368837a08b2a57d8b0671974af11756a53bb3a1a835b8ac5ecb3` | **SEMANTIC_CHANGE_REQUALIFIED** |
| `arrhenius_fracture/bulk_pt_live_run2d_composition_v1043.py` | `write_joint_veto_audit / post-rollback audit` | `61e34a2a61cce804d0490586993411ff1a8cb680d5e175e074280b7df5c510d6` | `61e34a2a61cce804d0490586993411ff1a8cb680d5e175e074280b7df5c510d6` | **IDENTICAL** |

Corrected-lineage bounded regression: **2 passed, 4 deselected**. No trajectory was run.

Lifecycle result: **QUALIFIED_FROM_EXISTING_TRANSACTIONAL_EVIDENCE**.
