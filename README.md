[![CI](https://github.com/gap-packages/GaloisGroups/actions/workflows/CI.yml/badge.svg)](https://github.com/gap-packages/GaloisGroups/actions/workflows/CI.yml)
[![Code Coverage](https://codecov.io/github/gap-packages/GaloisGroups/coverage.svg?branch=master&token=)](https://codecov.io/gh/gap-packages/GaloisGroups)


# The GAP 4 package `GaloisGroups'

This package implements the computation of Galois group.

For arithmetic computations, the package relies on
[PARI/GP](https://pari.math.u-bordeaux.fr/)  called via the PARIInterface
package.

## Requirements

GaloisGroups depends on other Galois packages

- GAPDoc
- Digraphs
- ferret
- alnuth
- TransGrp
- PARIInterface

as well as the development version of PARI/GP (see Installation below)

## Installation

For the moment you need to use the development version of pari
that implement function unavailable in PARI/GP 2.11 (e.g. `ZpXQX_liftroots`).
If you need the compiler to know about a special location of PARI/GP runs
configure with appropriate options

    $ ./configure [GAPPATH] [--with-pari=PARIPATH]


## Documentation

Full information and documentation can be found in the manual, available
as PDF `doc/manual.pdf` or as HTML `doc/chap0_mj.html`, or on the package
homepage at

  <https://gap-packages.github.io/GaloisGroups/>


## Bug reports and feature requests

Please submit bug reports and feature requests via our GitHub issue tracker:

  <https://github.com/gap-packages/GaloisGroups/issues>


## License

GaloisGroups is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

For details see the files COPYRIGHT.md and LICENSE.

## Testing

Tests are implemented in tst/. To run them do

    $ {GAP} tst/testall.g

## References

- C. Fieker, J. Klüners "Computation of Galois groups of rational polynomials" (2014)
- A. Hulpke "Galois groups through invariant relations" (1999)
- A. Hulpke "Finding intermediate subgroups" (2017)
