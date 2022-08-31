# sound_prop

This example contains a simple sound propagation simulator using `libjmm`. It consists of a few command-line programs which can be run to propagate sound and process the output.

## `mesh`

## `simulate`

## `render`

## `slice`

## Debugging

### viz.py

This script can be used to quickly visualize what's going with `simple_building.c`. While debugging, try doing the following:

1. Call `mesh3_dump_verts` and `mesh3_dump_cells` to write `verts.bin`
   and `cells.bin` to disk. Make sure that `viz.py` is pointing to
   these two files.
2. Likewise, call `eik3_dump_jet` to write `jet.bin` to disk.
3. Run `viz.py` (ideally from IPython in Emacs for more interactivity).
4. Call `mesh3_get_mean_edge_length` to get a good characteristic
   length scale for the mesh.
