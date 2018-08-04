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

The algorithm used for modeling propagation is [Longley-Rice]. Splash contains
[an implementation] that was initially ported from the C++ source but thereafter
modified to allow concurrent use. But the best thing about our implementation is
that it is entirely documented and brought back to understandable and meaningful
names, as well as containing references, descriptions, and explanations.

As such, I believe it is much more approachable than the C++ or FORTRAN versions
and that one does not need to refer back to the memos describing the algorithm.

While this tool and library is licensed, its Longley-Rice implementation and
associated documentation is released in the Public Domain.

[Longley-Rice]: https://en.wikipedia.org/wiki/Longley%E2%80%93Rice_model
[an implementation]: https://docs.rs/splash/*/splash/itm/index.html
