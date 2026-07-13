# Arrhenius FEM/CZM MPZ v9.1 three-class delta

This repository contains the verified update from the previously delivered
`Arrhenius_FEM_CZM_MPZ_v9_0` package to the three-class moving-process-zone
v9.1 tuning workflow.

The core payload adds or replaces **18 source, test, configuration, target, and
documentation files**. It does not delete any v9.0 files, and it omits only
regenerable smoke-output tables and figures.

## Obtain the repository

```bash
git clone https://github.com/ukaiiaku-maker/Arrhenius_FEM_CZM_MPZ.git
cd Arrhenius_FEM_CZM_MPZ
```

## Apply without overwriting v9.0

```bash
python apply_mpz_v9_1_delta.py \
  --base /path/to/Arrhenius_FEM_CZM_MPZ_v9_0 \
  --out /path/to/Arrhenius_FEM_CZM_MPZ_v9_1_three_class_tuning
```

If `--out` is omitted, the script creates a sibling folder named
`Arrhenius_FEM_CZM_MPZ_v9_1_three_class_tuning`.

Then run:

```bash
cd /path/to/Arrhenius_FEM_CZM_MPZ_v9_1_three_class_tuning
conda activate arrhenius-fem-czm
python -m pip install -e .
pytest -q tests/test_moving_process_zone.py tests/test_mpz_three_class_fit.py
STAGE=smoke bash run_mpz_three_class_tuning.sh
```

## Payload verification

```text
SHA-256: 8a49221c62faa9101e7fdc50968ddcb99aed3d74a24fd2f92885772a55126102
Files: 18
Base64 parts: 4
Decoded ZIP size: 36,234 bytes
```

The installer verifies the checksum, ZIP integrity, payload file count, and
path safety before writing any files. By default, it copies v9.0 into a new
version-specific folder and leaves the original package unchanged.
