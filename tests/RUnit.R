if (requireNamespace("RUnit", quietly=TRUE) && requireNamespace("dynafluxr", quietly=TRUE)) {
  testSuite <- defineTestSuite(
    name = "dynafluxr unit tests",
    dirs = system.file("unitTests", package = "dynafluxr"),
    testFuncRegexp = "^[Tt]est.+"
  )
  Sys.setenv("R_TESTS"="")
  tests <- runTestSuite(testSuite)

  printTextProtocol(tests)

  if (getErrors(tests)$nFail > 0) stop("RUnit test failure")
  if (getErrors(tests)$nErr > 0) stop("Errors in RUnit tests")
}
