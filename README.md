# libjmm

This is a C library implementing a variety of *jet marching methods*
(JMMs) for solving the eikonal equation in two and three dimensions.

## Build dependencies
1. Meson
   - Can be installed via pip or conda
   ``` sh
   pip install meson
   ```
2. Cython
   - Also available via pip/conda
   ```sh
   pip install cython
   ```
2. Other build pipeline tools:
   ``` sh
   sudo apt install cmake ninja-build build-essential pkg-config libgsl-dev
   ```
   Notes:
   - Requires GCC-10 or higher: https://askubuntu.com/questions/466651/how-do-i-use-the-latest-gcc-on-ubuntu/1163021#1163021
   - If you want to remove the warnings about an older version of cmake, make sure you have cmake version > 3.17.  If apt isn't finding the latest version, follow the instructions here to add the apt repository directly from cmake/kitware: https://askubuntu.com/questions/355565/how-do-i-install-the-latest-version-of-cmake-from-the-command-line

## Building

This library can be built using [Meson](https://mesonbuild.com/).

To build, run the following:

``` sh
git clone https://github.com/sampotter/jmm
cd jmm
meson setup --buildtype=release builddir
cd builddir
meson compile
```

## Dependencies

The library currently depends on a fork of [TetGen](https://github.com/sampotter/tetgen), which we maintain. This dependency will be picked up and built automatically by Meson.

## Documentation

Until this library is officially released, **there is no official documentation**. Instead, we recommend either taking a look at [the examples](#examples), or submitting an issue. Eventually, we will release tutorial-style documentation which gives an overview of the library, discusses some of the main entrypoints, and provides a roadmap.

In the meantime, individual functions in the library are documented, with the documentation comments preceding function definitions (i.e., the implementations are documented in the `.c` files). If you use something like [LSP Mode](https://emacs-lsp.github.io/lsp-mode/), it should be possible to look this documentation up on the fly. For this reason, it is unlikely that Doxygen-style documentation will be provided for this library.

## Examples

The `./examples` directory contains some examples of short programs which use this library to do different things. Each example contains its own `README.md` file explaining how to compile and run it. These examples can be used to quickly learn how to get going with this library.

## Development

Main development takes place in the `develop` branch.

## Tagged versions

Some important versions are tagged (you can find these under the
"branches" menu if you're browsing this repository using
GitHub). Currently, there are the following branches:

| Tag    | Description                                      |
|--------|--------------------------------------------------|
| v0.1.1 | First version with ~eik~ working for s \neq 1.   |
| v0.1.0 | First version with ~eik~ working for s \equiv 1. |
