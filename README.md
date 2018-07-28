# Splash

**[SPLAT!]** is a wonderful tool for RF analysis over terrain. Unfortunately,
it only works with a few terrain sources, supports no more than four emitters
in one analysis, and makes use of only a single core.

**Splash** is not a rewrite of SPLAT!. It is a different implementation of the
_idea_ of SPLAT!, but taking different routes and making different decisions.

## Differences

 - Splash only understands metres, and internally, it uses SI units.
 - Splash calculates local geocentric radii instead of assuming spherical Earth.
 - Splash does not output pretty pictures, but raster GIS data.
 - Splash uses geodesic distance calculations ([Karney 2013](https://doi.org/10.1007/s00190-012-0578-z)).

