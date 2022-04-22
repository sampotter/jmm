#!/usr/bin/env python

import numpy as np
import subprocess

from pathlib import Path

exe_path = Path('../../build/examples/3d_wedge/3d_wedge')

# problem parameters
ns = [0.25]
maxvols = [0.005]
rfacs = [0.2]
phips = [np.pi/4]
sps = [np.sqrt(2)]
widths = [4]
heights = [2]

for n, maxvol, rfac, phip, sp, width, height in zip(
        ns, maxvols, rfacs, phips, sps, widths, heights):
    # run the experiment
    args = [
        exe_path,
        f'-n{n}',
        f'--maxvol={maxvol}',
        f'--rfac={rfac}',
        f'--phip={phip}',
        f'--sp={sp}',
        f'--width={width}',
        f'--height={height}',
        f'--verbose'
    ]
    print(subprocess.list2cmdline(args))
    result = subprocess.run(args, capture_output=True, text=True)

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
                     'direct_jet.bin', 'direct_jet_gt.bin',
                     'direct_state.bin', 'direct_origin.bin',
                     'direct_par_l.bin', 'direct_par_b.bin',
                     'o_refl_jet.bin', 'o_refl_jet_gt.bin',
                     'o_refl_state.bin', 'o_refl_origin.bin',
                     'o_refl_par_l.bin', 'o_refl_par_b.bin',
                     'n_refl_jet.bin', 'n_refl_jet_gt.bin',
                     'n_refl_state.bin', 'n_refl_origin.bin',
                     'n_refl_par_l.bin', 'n_refl_par_b.bin']:
        path = Path(path_str)
        path.replace(out_path/path)
