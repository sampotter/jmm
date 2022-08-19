# `make_L_video.py`

This script sets up a point source problem in an "L-shaped room" test
problem.

Before running, make sure that `L_verts.bin` and `L_cells.bin` are in
this directory. These are serialized numpy arrays containing the
vertex and cell connectivity information for the tetrahedron mesh, and
be made by running:

```sh
$ ./setup.sh
```

in this directory.
