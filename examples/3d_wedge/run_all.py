#!/usr/bin/env python

import numpy as np
import subprocess

from pathlib import Path

exe_path = Path('../../build/examples/3d_wedge/3d_wedge')

# problem parameters
ns = [0.5]
maxvols = [0.01, 0.0333, 0.001]
rfacs = [0.2]
phips = [3*np.pi/4]
sps = [np.sqrt(2)]
widths = [4]
heights = [2]

for n, maxvol, rfac, phip, sp, width, height in zip(
        ns, maxvols, rfacs, phips, sps, widths, heights):
    print(f'n={n}, maxvol={maxvol}, rfac={rfac}, ' +
          f'phip={phip}, sp={sp}, width={width}, ' +
          f'height={height}')

    # run the experiment
    result = subprocess.run([exe_path,
                             f'-n{n}',
                             f'--maxvol={maxvol}',
                             f'--rfac={rfac}',
                             f'--phip={phip}',
                             f'--sp={sp}',
                             f'--width={width}',
                             f'--height={height}',
                             f'--verbose'],
                            capture_output=True,
                            text=True)

    # make new directory to store results
    out_path = Path(f'n{n}_a{maxvol}_rfac{rfac}_phip{phip}_sp{sp}_w{width}_h{height}')
    out_path.mkdir()

    # write output to results directory
    with open(out_path/'output.txt', 'w') as f:
        print(result.stdout, file=f)

    # move all output binary files to results directory
    for path_str in ['spec.txt',
                     'verts.bin', 'cells.bin',
                     'surface_faces.bin', 'surface_verts.bin',
                     'direct_jet.bin', 'direct_jet_gt.bin', 'direct_state.bin',
                     'o_refl_jet.bin', 'o_refl_jet_gt.bin', 'o_refl_state.bin',
                     'n_refl_jet.bin', 'n_refl_jet_gt.bin', 'n_refl_state.bin']:
        path = Path(path_str)
        path.replace(out_path/path)
