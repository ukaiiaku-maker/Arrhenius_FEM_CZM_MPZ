#!/usr/bin/env python3
from __future__ import annotations
import inspect, subprocess, sys
from pathlib import Path

def main():
 import arrhenius_fracture
 from arrhenius_fracture import sharp_front
 from arrhenius_fracture import mixed_mode_first_passage_v1 as mm
 print('MIXED_MODE_V1 project package:', Path(inspect.getfile(arrhenius_fracture)).resolve())
 print('MIXED_MODE_V1 sharp_front:', Path(inspect.getfile(sharp_front)).resolve())
 print('MIXED_MODE_V1 module:', Path(inspect.getfile(mm)).resolve())
 print('MIXED_MODE_V1 model:', mm.MODEL_ID)
 help_text=subprocess.run([sys.executable,'-m','arrhenius_fracture.sharp_front','--help'],text=True,capture_output=True)
 if help_text.returncode: raise SystemExit(help_text.stderr)
 required=['--crack-backend','--stop-after-first-fire','--emit-barrier-kind','--cleave-barrier-kind','--rJ-outer','--max-fronts']
 missing=[x for x in required if x not in help_text.stdout+help_text.stderr]
 if missing: raise SystemExit('active sharp_front is missing required options: '+', '.join(missing))
 cp=subprocess.run([sys.executable,'-m','unittest','tests.test_mixed_mode_first_passage_v1'])
 if cp.returncode: raise SystemExit(cp.returncode)
 print('MIXED_MODE_V1 verification OK')
if __name__=='__main__':main()
