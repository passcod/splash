# Splash

**[SPLAT!]** is a wonderful tool for RF analysis over terrain. Unfortunately,
it only works with a few terrain sources and only uses one core for processing.

**Splash** is not a rewrite of SPLAT!. It is a different implementation of the
_idea_ of SPLAT!, but taking different routes and making different decisions.

[SPLAT!]: http://www.qsl.net/kd2bd/splat.html

## Differences

 - Splash only understands metres, and internally, it uses SI units.
 - Splash does not output pretty pictures, but raster GIS data (GeoTIFF).
 - Splash does not use SDF input files, but raster GIS data (GeoTIFF).

## Longley-Rice

The algorithm used for modelling propagation is [Longley-Rice]. Splash contains
[an implementation] that was initially ported from the C++ source but thereafter
modified to allow concurrent use. But the best thing about our implementation is
that it is entirely documented and brought back to understandable and meaningful
names, as well as containing references, descriptions, and explanations.

As such, I believe it is much more approachable than the C++ or FORTRAN versions
and that one may not need to refer back to the memos describing the algorithm.

While this tool and library is licensed, its Longley-Rice implementation and
associated documentation is released in the Public Domain.

[Longley-Rice]: https://en.wikipedia.org/wiki/Longley%E2%80%93Rice_model
[an implementation]: https://docs.rs/splash/*/splash/itm/index.html

### Parallel performance expectations

Using the reimplementation and Rust's capabilities, and taking prior research
by [Song 2011] and [Musselman 2013] into account, Splash expects to be able to
perform at least as fast given intervening years' compute advances.

Notably, **Song 2011** describes obtaining a 150k point "HD" splat propagation
model in less than a second on an nVidia Tesla GPU. Current consumer-grade
Intel GPUs are several times more powerful, and gaming-grade GPUs may be orders
of magnitude faster. Even an Intel i7-core CPU with symmetric multiprocessing
may bring comparable results. Thus, the expectation with Splash is to obtain, at
minimum:

 - with GPU: < 500 milliseconds for a full render, or
 - CPU only: < 20 seconds for a full render.

[Song 2011]: https://ieeexplore.ieee.org/document/5680900/
[Musselman 2013]: https://github.com/amusselm/Parallel-LRP
