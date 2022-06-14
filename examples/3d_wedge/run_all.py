#!/usr/bin/env python

import itertools as it
import numpy as np
import subprocess
import time

from pathlib import Path

exe_path = Path('../../build/examples/3d_wedge/3d_wedge')

# problem parameters
ns = [1.75]
p0, p1 = 2, 4
maxvols = np.logspace(-p0, -p1, 2*(p1 - p0) + 1)
rfacs = [0.3]
phips = [np.pi/4]
sps = [np.sqrt(2)]
widths = [4]
heights = [2]

for n, maxvol, rfac, phip, sp, width, height in it.product(
        ns, maxvols, rfacs, phips, sps, widths, heights):

    # make new directory to store results
    out_path = Path(
        f'n{n}_a{maxvol}_rfac{rfac}_phip{phip}_sp{sp}_w{width}_h{height}')

    # if the experiment already exists, skip it and warn
    if out_path.exists():
        print(f'experiment already exists!')
        continue
    out_path.mkdir()

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
    print(f'Running command: "{subprocess.list2cmdline(args)}"')

    t0 = time.time()
    result = subprocess.run(args, capture_output=True, text=True)
    print(f'- elapsed time: {time.time() - t0:.1f}s')

    # write output to results directory
    with open(out_path/'output.txt', 'w') as f:
        print(result.stdout, file=f)

    # move all output binary files to results directory
    for path_str in ['spec.txt',
                     'verts.bin', 'cells.bin',
                     'surface_faces.bin', 'surface_verts.bin',
                     'direct_jet.bin', 'direct_jet_gt.bin', 'direct_hess.bin',
                     'direct_state.bin', 'direct_origin.bin',
                     'direct_par_l.bin', 'direct_par_b.bin',
                     'direct_accepted.bin',
                     'direct_t_in.bin', 'direct_t_out.bin',
                     'o_refl_jet.bin', 'o_refl_jet_gt.bin', 'o_refl_hess.bin',
                     'o_refl_state.bin', 'o_refl_origin.bin',
                     'o_refl_par_l.bin', 'o_refl_par_b.bin',
                     'o_refl_accepted.bin',
                     'o_refl_t_in.bin', 'o_refl_t_out.bin',
                     'n_refl_jet.bin', 'n_refl_jet_gt.bin', 'n_refl_hess.bin',
                     'n_refl_state.bin', 'n_refl_origin.bin',
                     'n_refl_par_l.bin', 'n_refl_par_b.bin',
                     'n_refl_accepted.bin',
                     'n_refl_t_in.bin', 'n_refl_t_out.bin',
                     'img_grid.txt',
                     'slice_direct_A.bin', 'slice_direct_T.bin',
                     'slice_direct_origin.bin', 'slice_direct_E_T.bin',
                     'slice_o_refl_A.bin', 'slice_o_refl_T.bin',
                     'slice_o_refl_origin.bin', 'slice_o_refl_E_T.bin',
                     'slice_n_refl_A.bin', 'slice_n_refl_T.bin',
                     'slice_n_refl_origin.bin', 'slice_n_refl_E_T.bin']:
        path = Path(path_str)
        path.replace(out_path/path)
