This file describes changes in the GaloisGroups package.

## Unreleased

- Replace the `GaloisGroup` attribute by `Galois`, which computes Galois groups
  of polynomials using PARI through a kernel extension
- Add `GaloisDescentTable`, `PrintGaloisDescentTable` and precomputed descent
  tables `GaloisDescentTables` up to degree 17
- Add test polynomials `GaloisTestPolynomials` up to degree 17
- Speed up `AllMonomials` and the computation of flat monomials
- Fix the generation of resolvents and a mix-up of G and H in `Galois`
- Require GAP >= 4.11 and the PARIInterface and TransGrp packages

## 0.1 (2018-01-31)
