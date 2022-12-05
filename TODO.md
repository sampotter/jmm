# Todo

## Misc

* [ ] Inline PLY describing building test problem
  - don't need to load file
  - call TetGen on PLY directly
* [ ] Start working on documentation

## Refactoring `eik3`

* [ ] Rename `eik2` to `eik2f` (or something)
* [ ] Rename `eik3` to `eik31`
* [ ] Move `transport`-related stuff into a separate module
* [ ] Clean up `eik3`'s API
* [ ] Combine `jet3` and `hess`
* [ ] Simplify logic in `update`
* [ ] Track `t_in` and `t_out` separately... move into separate module
      along with transport stuff

## Improving the `eik3` algorithm

* [ ] Get rid of diffraction `ftype`

## Results

* [ ] Ray-traced video of direct arrival wavefront propagating
      throughout building with opacity selected by amplitude and
      colored by error
  - coloring by error is tough but doable... one way to go would be to
    do it automatically would be to use Richardson extrapolation...
* [ ] Fully worked out wedge problem with errors for all quantities

# Some notes and observations

As much as possible, we want to get rid of our extra stuff related to
updating from edges. All the `bde` functions... We don't want to set
BCs here, since we don't want to do "diffraction" updates.

What's the right way to think about updating? At each step, we want to
minimize over the VALID front. Probably the best way to do this is to
just do so explicitly. But considering the triangle fan around the
newly VALID node makes things a bit simpler.

Pretty sure we want to completely get rid of t_out.

Probably would be good to get rid of `do_diff_edge_updates_and_adjust`.

# Lessons learned

* Explicit error handling from the get-go

# JCP revision

- [X] factor out caches
- [ ] introduce sfunc
