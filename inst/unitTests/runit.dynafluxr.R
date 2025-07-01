#library(RUnit)
norm2=function(x) sum(x*x)
test.cli.ex=function() {
  ddir=system.file("dataglyco", package="dynafluxr")
  meas=file.path(ddir, "data.tsv")
  sto=file.path(ddir, "network.txt")
  res=dynafluxr::cli(c("-m", meas, "-s", sto, "--skip", "24", "-o", ""))
  residm=res$rsp(res$mf$Time)
  #print(norm2(residm))
  RUnit::checkTrue(norm2(residm) < 1.e-4)
}
