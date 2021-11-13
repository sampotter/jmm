# `visualize_cutset.py`

This script solves a point source problem for a constant speed of
sound in an L-shaped room. Afterwards, it extracts the shadow cutset,
plotting the edges in the cutset, the approximated shadow surface, and
the nodes in the tetrahedron mesh. Green nodes are `valid`, purple
nodes are `shadow`, and the highlighted red nodes is the original
point source. The point where each cutset intersects the shadow
surface is indicated with a white sphere. Note that the solver stops
short when it reaches the shadow zone to avoid computing unnecessary
information there.

After running, a quick "preview" screenshot will be written to
`<outpath>/cutset.png`, and a
[vtk.js](https://kitware.github.io/vtk-js/index.html) scene will be
written to `<outpath>/cutset.vtkjs`. To preview the scene, visit
`viewer.pyvista.org` and drag and drop the vtkjs file there.

Before running, make sure that `L_verts.bin`, `L_cells.bin`, and
`L.obj` are in this directory by running:

```sh
$ ./setup.sh
```

in this directory.

If you run `./run_for_all_point_sources.sh`, it will iterate through
each vertex in the mesh, place a point source there, run
`visualize_cutset.py`, and store the output in `srcs/<point source
index>`.
