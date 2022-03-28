library(RUnit)
norm2=function(x) sum(x*x)
test.cli.ex=function() {
  ddir=system.file("data", package="dynafluxr")
  meas=file.path(ddir, "data_teusink.tsv")
  sto=file.path(ddir, "network_teusink.txt")
  res=dynafluxr::cli(c("-m", meas, "-s", sto, "--skip", "10"))
  tp=res$mf$Time
  np=length(tp)
  residm=res$rsp(tp)
  checkTrue(norm2(residm) < 1.e-21)
}
