libpredict
==========

A satellite orbit prediction library.


Building
--------

We recommend using out-of-source builds.

```
cd $SOURCEDIR
mkdir build
cd build
cmake ..
make
```


Installation
------------

```
make install
```

The install location is defined by `CMAKE_INSTALL_PREFIX`, which
defaults to `/usr/local`. To relocate the whole installation (to make
usr/local etc. inside another directory, e.g., if you want to make a
tarball or package it afterwards), use `make DESTDIR=/foo/bar install`.

Linking
-------

The library comes with pkg-config information, so the include and
library paths and flags can be found using the `pkg-config` command or
by using the `PKG_CHECK_MODULES` autotools macro or CMake command.

Quickstart
----------

We recommend to investigate the examples under `examples/` and the API documentation in `include/predict/predict.h`.

A condensed version is that
```
predict_orbital_elements_t *orbital_elements = orbital_elements = predict_parse_tle(tle_line_1, tle_line_2);
```
parses a TLE in the form of two char arrays to orbital elements representing a satellite, and that
```
predict_observer_t *observer = predict_create_observer(name, latitude_radians, longitude_radians, altitude_meters);
```
defines a QTH for observation. For prediction,
```
struct predict_position orbit;
predict_orbit(orbital_elements, &orbit, prediction_time);
```
can be used to calculate properties that are independent of an observer (longitude, latitude, ...), while
```
struct predict_observation observation;
predict_observe_orbit(observer, &orbit, &observation);
```
will convert to properties that will be relative to our observer (azimuth, elevation, ...).

License
-------

 Copyright 1991-2006 John A. Magliacane (KD2BD)
 
 Copyright 2013- Akademisk radioklubb (LA1K)
 
 Copyright 2013-2015 Knut Magnus Kvamtrø (LA3DPA)

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
