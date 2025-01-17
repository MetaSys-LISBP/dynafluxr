## v0.20.0 2025-01-17
- all NA species are now removed from stoichiometrique matrix
- curvature of L-curve is used for ref knot number
- sf is now applied to var_ref too.
- point character is now "." instead of "o"
- line width is set to 1.5 instead of 1
- stricter separation of monotonicity for preliminary fit and sto fit
- NA in fits gives an error now
- added stoinv to res
- stricter control of negative coefficients in B-splines

## v0.19.0 2024-12-06
- sf in file
- sf for all => error
- total chi2
- manual SD
- res$internal_knot_ref
- example of undefined rates

## v0.18.0 2024-12-03
- added minimal norm for under-determined stoichiometric matrix

## v0.17.0 2024-10-28
- added --sf to ILS
- passed to regular knot grid for bspline

## v0.16.0 2024-10-21
- converted --wsd to use of gmresls in DLS
