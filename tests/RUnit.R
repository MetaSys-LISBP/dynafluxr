if (requireNamespace("RUnit", quietly=TRUE) && requireNamespace("dynafluxr", quietly=TRUE)) {
  testSuite <- RUnit::defineTestSuite(
    name = "dynafluxr unit tests",
    dirs = system.file("unitTests", package = "dynafluxr"),
    testFuncRegexp = "^[Tt]est.+"
  )
  Sys.setenv("R_TESTS"="")
  tests <- RUnit::runTestSuite(testSuite)

  RUnit::printTextProtocol(tests)

  if (RUnit::getErrors(tests)$nFail > 0) stop("RUnit test failure")
  if (RUnit::getErrors(tests)$nErr > 0) stop("Errors in RUnit tests")
}
