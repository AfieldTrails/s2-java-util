
# S2 Geometry utilities

This library is built on [Google's S2 geometry library](https://github.com/google/s2-geometry-library-java), providing some optimizations (`S2CellIds`) and utilities (everything else).

There are several excellent blog posts describing what Google's library provides:
  * [Google’s S2, geometry on the sphere, cells and Hilbert curve](http://blog.christianperone.com/2015/08/googles-s2-geometry-on-the-sphere-cells-and-hilbert-curve/) by Christian S. Perone
  * [S2 cells and space-filling curves: Keys to building better digital map tools for cities](https://medium.com/sidewalk-talk/s2-cells-and-space-filling-curves-keys-to-building-better-digital-map-tools-for-cities-a312aa5e2f59) by Sven Kreiss

# `S2CellIds`

This version of `S2CellId` has been optimized for use Android to reduce the number of allocations required. It can be used to replace `S2CellId` objects with `long`.

# `S2Distance`

Has some utility methods for working with distances in the S2 library. Most notably, it contains a port of the Android Location distance method that computes the distance on the WGS84 spheroid. This allows distance calculations to more closely match between representations.
