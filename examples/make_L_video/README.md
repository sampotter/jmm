# `make_L_video.py`

This script sets up a point source problem in an "L-shaped room" test problem. Specifically,  $\Omega = ([0, 2] \times [0, 2] \times [0, 1]) \backslash ([1, 2] \times [1, 2] \times [0, 1])$.

Before running, make sure that ~L_verts.bin~ and ~L_cells.bin~ are in this directory. These are serialized numpy arrays containing the vertex and cell connectivity information for the tetrahedron mesh, and be made by running:

```sh
$ ./setup.sh
```

in this directory.
