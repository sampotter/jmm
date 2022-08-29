[![Build Status](https://app.travis-ci.com/sampotter/jmm.svg?branch=master)](https://app.travis-ci.com/sampotter/jmm)

# libjmm

This is a C library implementing a variety of *jet marching methods*
(JMMs) for solving the eikonal equation in two and three dimensions.

## Building

This library can be built in the usual way with [Meson](https://mesonbuild.com/) or [CMake](https://cmake.org/).

Initially only CMake was supported, but we are switching to use Meson as the primary build system. The only `CMakeLists.txt` files will be kept for now.

The library currently depends on a fork of [TetGen](https://github.com/sampotter/tetgen), which we maintain. If you're using Meson to build libjmm, this dependency will be picked up automatically, and there is nothing extra that needs to be done. For CMake, you must clone a copy of this repository to `./subprojects/tetgen`.

## Documentation

Until this library is officially released, there will not be any documentation. Instead, we recommend either taking a look at [the examples](#examples), or submitting an issue. Eventually, we will release tutorial-style documentation which gives an overview of the library, discusses some of the main entrypoints, and provides a roadmap.

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
