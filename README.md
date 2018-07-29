# Splash

**[SPLAT!]** is a wonderful tool for RF analysis over terrain. Unfortunately,
it only works with a few terrain sources and only uses one core for processing.

**Splash** is not a rewrite of SPLAT!. It is a different implementation of the
_idea_ of SPLAT!, but taking different routes and making different decisions.

## Differences

 - Splash only understands metres, and internally, it uses SI units.
 - Splash calculates local geocentric radii instead of assuming spherical Earth.
 - Splash uses Vicenty distance calculations rather than Haversine most times.
 - Splash does not output pretty pictures, but raster GIS data (GeoTIFF).
 - Splash does not use SDF input files, but raster GIS data (GeoTIFF).
