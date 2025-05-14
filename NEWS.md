## v0.26.0 2025-05-14
- added legal information to gui()
- added error for mono in DLS mode
- gui(): in case of error, results are canceled

## v0.25.0 2025-04-16
- added 'date' and a short 'call' for R and Shell (without default params)
- added scroll to top at the end of calculations
- launch browser even in script mode

## v0.24.0 2025-03-14
- added gui()

## v0.23.0 2025-03-05
- added option `regular_grid` (default: TRUE)
- added date and opt to res

## v0.22.1 2025-02-07
- added isp() SD for --wsd too

## v0.22.0 2025-02-06
- added SD for isp() (no --wsd yet)

## v0.21.0 2025-01-22
- added constraints by value at given time points for ILS
- added option "--pch", default to "."
- added nm_lapply()
- added res$mffull
- curvature of L-curve is measured on log-values
- fixed error message for non measured metabs but asked with scaling factors

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
