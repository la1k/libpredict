libpredict
==========

A satellite tracking library.

Building
--------

We recommend using out-of-source builds, e.g., by creating a `build` directory in the source tree:

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

The install location is defined by `CMAKE_INSTALL_PREFIX`, which defaults to `/usr/local`. To relocate the whole installation (to make usr/local etc. inside another directory, e.g., if you want to make a tarball or package it afterwards), use `make DESTDIR=/foo/bar install`.
