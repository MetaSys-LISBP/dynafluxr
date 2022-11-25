# dynafluxr

Reaction rate dynamics can be retrieved from time-series measurements
  of chemical species and/or their isotopes. User has to provide
  corresponding stoechiometric
  matrix but not a regulation model (Michaelis-Menten or similar).
  Instead of solving an ODE system describing the evolution of concentrations,
  we use B-splines to catch the concentration and rate dynamics then solve a least square problem
  on their coefficients with optional constraints. Constraints can be set
  on initial values of concentration and/or monotonicity of B-spline coefficients.
  Positivity of these coefficients is always enforced.
  The package 'dynafluxr' can be used as a library providing appropriate functions but also as
  a stand alone application with command line interface.
