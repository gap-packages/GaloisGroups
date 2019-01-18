[![Build Status](https://travis-ci.org/gap-packages/GaloisGroups.svg?branch=master)](https://travis-ci.org/gap-packages/GaloisGroups)
[![Code Coverage](https://codecov.io/github/gap-packages/GaloisGroups/coverage.svg?branch=master&token=)](https://codecov.io/gh/gap-packages/GaloisGroups)


# The GAP 4 package `GaloisGroups'

An implementation of Fieker and Klüners algorithm to compute galois groups
as described in https://arxiv.org/pdf/1211.3588.pdf

For arithmetic computations it relies on PARI/GP called via the
PARIInterface package.

## Documentation

Full information and documentation can be found in the manual, available
as PDF `doc/manual.pdf` or as HTML `htm/chapters.htm`, or on the package
homepage at

  <http://gap-packages.github.io/GaloisGroups/>


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

## TODO

- given a transitive subgroup G, we need to iterate over all subgroups
  in between Stabilizer(G,1) and G.

- TuplesSetsTuplesSets action (corresponding to permutation action
  on multivariate polynomials)

