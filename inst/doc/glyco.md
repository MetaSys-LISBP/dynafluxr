---
title: "Metabolic Flux Dynamics in Glycolysis with `dynafluxr`"
output: rmarkdown::html_vignette
author:
  - "Serguei Sokol^[TBI, Université de Toulouse, CNRS, INRAE, INSA, Toulouse, France] ^[MetaboHUB, National Infrastructure of Metabolomics and Fluxomics, Toulouse 31077, France]"
date: 2022-11-28
vignette: >
  %\VignetteIndexEntry{glyco}
  %\VignetteEngine{knitr::knitr}
  %\VignetteEncoding{UTF-8}
---




``` r
library(dynafluxr)
```
## Introduction

This vignette shows how `dynafluxr` R package can be used to unravel metabolic flux dynamics from metabolite time-series measurements and a stoichiometric model. No regulation model (like Michaelis-Menten) is required for this task. Our method is heavily relying on B-splines as implemented in `bspline` package. We show how to prepare measurement data, edit stoichiometric reactions, use various package options and interpret/exploit the results. The data used for this vignette were obtained by NMR measurements on MetaToul-FluxoMet platform, Toulouse Biotechnology Institute (TBI), Toulouse, France. Credits for making these data available go to Pauline Rouane, Cyril Charlier, Guy Lippens. Pierre Millard has shared stoichiometric model. All cited persons are from TBI.

## Data preparations

Data must be stored in a `.tsv` (Tabulation Separated Values) file (`.tsv` extension is not mandatory, it can be `.txt`, `.csv`, whatever, but the field separator must be the tabulation character), one column per chemical specie. The first row describes column names. Time points, at which measurements are done, must be stored in a column whose name starts with 'Time'. No specie name can start with 'Time'. Decimal separator must be a point `.` character. A head of data file used as an example is looking like


``` r
fmeas=system.file("dataglyco/data_PRO.tsv", package="dynafluxr")
cat(head(readLines(fmeas)), sep="\n")
#> Time	ACE	 ALA 	 FOR 	 PYR 	 ETOH 	 GLC 	 LAC 
#> 5	0.686079159909981	0.352668904657515	0.271293771262705	0.282469139463095	0.248351940239224	5.16630814568026	0.0712406341958987
#> 6.75	0.687672444646716	0.361318637275728	0.27391401585408	0.329132846916247	0.251620459437311	5.3910661637285	0.0716638244637927
#> 8.5	0.690616757561834	0.375205904823668	0.268544625230093	0.383292002434794	0.259032951308249	4.82270369767996	0.0762050140643337
#> 10.25	0.688103794715911	0.382373449436483	0.271928064222974	0.426170418853683	0.260870660454035	4.63965858112673	0.076885446185115
#> 12	0.69295232688529	0.403301241274558	0.266181057154522	0.489500439660487	0.271136289284411	4.59847699564914	0.0848684307472358
```

The raw data may be difficult to read because of tab misalignment and row wrapping. Here are parsed data:

``` r
mf=read.delim(fmeas, comment.char="#")
print(head(mf))
#>    Time       ACE       ALA       FOR       PYR      ETOH      GLC        LAC
#> 1  5.00 0.6860792 0.3526689 0.2712938 0.2824691 0.2483519 5.166308 0.07124063
#> 2  6.75 0.6876724 0.3613186 0.2739140 0.3291328 0.2516205 5.391066 0.07166382
#> 3  8.50 0.6906168 0.3752059 0.2685446 0.3832920 0.2590330 4.822704 0.07620501
#> 4 10.25 0.6881038 0.3823734 0.2719281 0.4261704 0.2608707 4.639659 0.07688545
#> 5 12.00 0.6929523 0.4033012 0.2661811 0.4895004 0.2711363 4.598477 0.08486843
#> 6 13.75 0.6963948 0.4200223 0.2626969 0.5508043 0.2774260 4.538150 0.09194118
```
## Stoichiometric model preparation

The reactions are written as in following example:

``` r
fsto=system.file("dataglyco/network_PRO.txt", package="dynafluxr")
cat(head(readLines(fsto)), sep="\n")
#> hxk	GLC + ATP -> G6P + ADP
#> pgi	G6P -> F6P
#> pfk	F6P + ATP -> FBP + ADP
#> fba	FBP -> DHAP + GAP
#> tpi	DHAP -> GAP
#> gapdh	GAP + NAD -> BPG + NADH
```
The first field separated by tabulation is a reaction name then reaction itself where species are separated by " + " sign and possibly preceded by a stoichiometric coefficient with " * " symbol.

We start by gathering example files in a freshly created working directory, say `glyco/`:

``` r
dir.create("glyco")
#> Warning in dir.create("glyco"): 'glyco' existe déjà
file.copy(c(fmeas, fsto), "glyco/")
#> [1] FALSE FALSE
```
## Data exploring

Let fit available data with B-splines of order 4 and with 5 internal knots (default in `dynafluxr`), then plot them to have an idea how data looks like:

``` r
fmeas="glyco/data_PRO.tsv"
mf=read.delim(fmeas, comment.char="#")
nm=colnames(mf)[-1L]
msp=fitsmbsp(mf$Time, mf[, nm], n=4, nki=5)
matplot(mf$Time, msp(mf$Time), t="l", ylab="Concentration [mM]", xlab="Time [s]")
matpoints(mf$Time, mf[, nm], pch="o", cex=0.5)
legend("topright", legend=nm, lty=1:5, col=1:6, cex=0.75)
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png)

We see that glucose (GLC) is going too close to 0, maybe even getting negative. A closer look confirms this (negative values are plotted in blue):

``` r
ineg=msp(mf$Time, "GLC") < 0
matplot(mf$Time[!ineg], msp(mf$Time, "GLC")[!ineg], t="l", ylab="C", xlab="Time", main="GLC")
abline(h=0)
matplot(mf$Time[ineg], msp(mf$Time, "GLC")[ineg], t="l", col="blue", lwd=2, add=TRUE)
matpoints(mf$Time, mf[,"GLC"], pch="o", cex=0.5)
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)

As concentrations cannot get negative values by definition, we constraint B-splines to non negative values with option `positive=1`

``` r
msp=fitsmbsp(mf$Time, mf[,-1], n=4, nki=5, positive=1)
ineg=msp(mf$Time, "GLC") < 0
print(sum(ineg))
#> [1] 0
matplot(mf$Time, msp(mf$Time, "GLC"), t="l", ylab="C", xlab="Time", main="GLC")
abline(h=0)
matpoints(mf$Time, mf[,"GLC"], pch="o", cex=0.5)
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png)

Negative values have disappeared but at the end of time interval, around 300 s, GLC curve is increasing which has no physical meaning. Glucose is only consumed in this experimental setup. Fortunately, we can constraint B-splines to be monotonously decreasing for GLC:


``` r
vmono=rep(0, length(nm))
vmono[match("GLC",nm)]=-1
msp=fitsmbsp(mf$Time, mf[,-1], n=4, nki=5, positive = 1, monotone = vmono)
matplot(mf$Time, msp(mf$Time, "GLC"), t="l", ylab="C", xlab="Time", main="GLC")
matpoints(mf$Time, mf[,"GLC"], pch="o", cex=0.5)
```

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png)

Looks good now. In `dynafluxr`, we have no option to constraint B-splines to non negative values, because it is done automatically. On the other hand, we do need to indicate which specie is supposed to be monotonously decreasing or increasing. It can be be done with `--mono MONO` option where `MONO` is a file name appropriately formatted to indicate which metabolite is increasing and which is decreasing. It can also be done with options `--incresing` and `--decreasing` with corresponding coma separated lists of specie names. It is this last option that we will use with `dynafluxr` but for now, let check some other metabolites, e.g. pyruvate (PYR) to be sure that the fit looks good too:

``` r
matplot(mf$Time, msp(mf$Time, "PYR"), t="l", ylab="C", xlab="Time", main="PYR")
matpoints(mf$Time, mf[,"PYR"], pch="o", cex=0.5)
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png)

Looks OK, so we are set to move to the next step.

## First estimate of fluxes with `dynafluxr`

We try `dynafluxr::cli()` (Command Line Interface) function with minimal options, just indicating measurement data, stoichiometric model and, as seen before, that GLC must be monotonously decreasing (`--decreasing` is abbreviated to `--decr`):


``` r
res=cli(c("-m", "glyco/data_PRO.tsv", "-s", "glyco/network_PRO.txt", "--decr", "GLC"))
```

The same command could be run from a shell (not R) environment. In this case it would look like:

```
$ Rscript -e "dynafluxr::cli()" -m glyco/data_PRO.tsv -s glyco/network_PRO.txt --decr GLC
```

The full list of available options can be viewed with `dynafluxr::cli("-h")` command. The result of `cli()` function (in case of no error in run-time), is a series of data and pdf files written in a directory `glyco/data_PRO` which is the stump of measurement data file:


``` r
list.files("glyco/data_PRO/")
#>  [1] "env.RData"   "flux.pdf"    "flux.tsv"    "ispecie.pdf" "ispecie.tsv"
#>  [6] "rate.pdf"    "rate.tsv"    "Readme.md"   "resid.pdf"   "specie.pdf" 
#> [11] "specie.tsv"  "stats.tsv"
```

The file `Readme.md` describes the content of the result directory. These files can be explored by system tools like pdf-viewers and spreadsheet software but for needs of this vignette we will reproduce some results by R means. As we won't need result files in this demo page, we'll cancel their writing with empty option `-o ""` which normally indicates result directory/archive name if it is different from the default one:

``` r
res=cli(c("-m", "glyco/data_PRO.tsv", "-s", "glyco/network_PRO.txt", "--decr", "GLC", "-o", ""))
```

Let now glance on fit quality, i.e. results of $\chi^2$-test:

``` r
print(res$chi2tab)
#>            rss      var_ref  df        chi2          pval
#> ACE   9.372441 3.910221e-06 180 2396908.128  0.000000e+00
#> ALA   2.227892 6.490101e-06 180  343275.469  0.000000e+00
#> FOR   2.100474 1.201390e-05 180  174836.897  0.000000e+00
#> PYR   2.230509 2.025976e-05 180  110095.514  0.000000e+00
#> ETOH  4.247570 1.392246e-05 180  305087.563  0.000000e+00
#> GLC  16.385392 8.677776e-03 180    1888.202 3.837689e-282
#> LAC   4.955555 2.384174e-05 180  207852.085  0.000000e+00
```

We can see that all $p$-values are 0 which means a very poor fit quality. Let see what is going on by plotting integrated species with measurements:


``` r
nm=colnames(res$mf)[-1]
matplot(res$tpp, res$isp(res$tpp, nm), t="l", xlab="Time", ylab="C")
matpoints(res$tp, res$mf[,nm], pch="o", cex=0.5)
legend("topright", legend=nm, lty=1:5, col=1:6, cex=0.75)
```

![plot of chunk unnamed-chunk-15](figure/unnamed-chunk-15-1.png)

Indeed, for all species, integrated curves are not fitting well the corresponding data. We take a closer look on pyruvate (PYR) which is buried in the full graph:

``` r
plot(res$tpp, res$isp(res$tpp, "PYR"), t="l", xlab="Time", ylab="C", ylim=range(res$mf$PYR))
points(res$tp, res$mf$PYR, pch="o", cex=0.5)
```

![plot of chunk unnamed-chunk-16](figure/unnamed-chunk-16-1.png)
We will use PYR example for tracking the fit quality as fit options change.

## Getting fit better

Let see the integral residual of metabolites which are not measured and thus are supposed to be close to 0:


``` r
nmall=rownames(res$sto)
nm=colnames(res$mf)[-1L]
nm0=setdiff(nmall, nm)
matplot(res$tpp, res$risp(res$tpp, nm0), t="l", ylab="M-\u222bSv dt", xlab="Time")
legend("topright", legend=nm0, lty=1:5, col=1:6, cex=0.75)
```

![plot of chunk unnamed-chunk-17](figure/unnamed-chunk-17-1.png)

We can see that the most diverging metabolites are fructose bi-phosphate (FBP) and fructose-6-phosphate (F6P). We can try to declare them as not available (NA) thus they won't degrade least-squares objective function with its diverging from 0 residuals:

``` r
res=cli(c("-m", "glyco/data_PRO.tsv", "-s", "glyco/network_PRO.txt", "--decr", "GLC", "-o", "", "--lna", "FBP,F6P"))
nmall=rownames(res$sto)
nm=colnames(res$mf)[-1L]
nm0=setdiff(nmall, nm)
layout(t(1:2))
matplot(res$tpp, res$risp(res$tpp, nm0), t="l", ylab="M-\u222bSv dt", xlab="Time")
legend("topright", legend=nm0, lty=1:5, col=1:6, cex=0.75)
plot(res$tpp, res$isp(res$tpp, "PYR"), t="l", xlab="Time", ylab="C", ylim=range(res$mf$PYR), main="PYR")
points(res$tp, res$mf$PYR, pch="o", cex=0.5)
```

![plot of chunk unnamed-chunk-18](figure/unnamed-chunk-18-1.png)
There is some progress in fit quality but PYR is still not well fitted. We try to add next diverging metabolites to NA: BPG (a little bit buried in the plot) and ACCOA.

``` r
res=cli(c("-m", "glyco/data_PRO.tsv", "-s", "glyco/network_PRO.txt", "--decr", "GLC", "-o", "", "--lna", "FBP,F6P,BPG,ACCOA"))
nmall=rownames(res$sto)
nm=colnames(res$mf)[-1L]
nm0=setdiff(nmall, nm)
layout(t(1:2))
matplot(res$tpp, res$risp(res$tpp, nm0), t="l", ylab="M-\u222bSv dt", xlab="Time")
legend("topright", legend=nm0, lty=1:5, col=1:6, cex=0.75)
plot(res$tpp, res$isp(res$tpp, "PYR"), t="l", xlab="Time", ylab="C", ylim=range(res$mf$PYR), main="PYR")
points(res$tp, res$mf$PYR, pch="o", cex=0.5)
```

![plot of chunk unnamed-chunk-19](figure/unnamed-chunk-19-1.png)

Now, residuals for "0" metabolites are all very close to 0 which means that they are not degrading least squares anymore and PYR is very well fitted. A good overall fit is confirmed by $\chi^2$-test:

``` r
print(res$chi2tab)
#>               rss      var_ref  df chi2      pval
#> ACE  0.0007429421 3.910221e-06 180  190 0.2902766
#> ALA  0.0012331192 6.490101e-06 180  190 0.2902766
#> FOR  0.0022826417 1.201390e-05 180  190 0.2902766
#> PYR  0.0038493543 2.025976e-05 180  190 0.2902766
#> ETOH 0.0026452678 1.392246e-05 180  190 0.2902766
#> GLC  1.6487773464 8.677776e-03 180  190 0.2902766
#> LAC  0.0045299301 2.384174e-05 180  190 0.2902766
```

as all $p$-values are well above 0.05 usual threshold.
We are ready to see the main result: estimated reaction rates.

``` r
matplot(res$tpp, res$vsp(res$tpp), t="l", ylab="Rate [1/s]", xlab="Time")
legend("topright", colnames(bsppar(res$vsp)$qw), lty=1:5, col=1:6, cex=0.75)
```

![plot of chunk unnamed-chunk-21](figure/unnamed-chunk-21-1.png)
It can be observed some oscillations in reaction rates and we must check if they don't come from a possible overfitting.

## Exploring problem of overfitting

By default, the number of internal knots used for fitting is 5. Let try few values under and above (1 to 6) and see if the measurements are still well fitted. We will continue to use PYR for this purpose:

``` r
layout(matrix(1:6, nrow=2, byrow=TRUE))
resk=lapply(1:6, function(k) {
  r=cli(c("-m", "glyco/data_PRO.tsv", "-s", "glyco/network_PRO.txt", "--decr", "GLC", "-o", "", "--lna", "FBP,F6P,BPG,ACCOA", "-k", k))
  plot(r$tpp, r$isp(r$tpp, "PYR"), main=paste0("PYR k=",k,sep=""), ylab="C", xlab="Time", t="l")
  points(r$tp, r$mf$PYR, pch="o", cex=0.5)
  r
})
```

![plot of chunk unnamed-chunk-22](figure/unnamed-chunk-22-1.png)
Starting from value `k=3`, the fits look satisfactory. Let see what says $\chi^2$ test:

``` r
for (k in 1:6) {
  cat("k=",k,"\n")
  print(resk[[k]]$chi2tab)
}
#> k= 1 
#>               rss      var_ref  df      chi2          pval
#> ACE  0.0009204292 3.910221e-06 184  235.3906  6.251096e-03
#> ALA  0.0013841113 6.490101e-06 184  213.2650  6.867265e-02
#> FOR  0.0023114324 1.201390e-05 184  192.3964  3.206649e-01
#> PYR  0.0992589090 2.025976e-05 184 4899.3133  0.000000e+00
#> ETOH 0.0027174157 1.392246e-05 184  195.1821  2.722272e-01
#> GLC  2.3669580284 8.677776e-03 184  272.7609  2.307959e-05
#> LAC  0.0271184672 2.384174e-05 184 1137.4367 4.456405e-137
#> k= 2 
#>               rss      var_ref  df      chi2          pval
#> ACE  0.0008863323 3.910221e-06 183  226.6706  1.550075e-02
#> ALA  0.0013736301 6.490101e-06 183  211.6500  7.210724e-02
#> FOR  0.0023144314 1.201390e-05 183  192.6461  2.979457e-01
#> PYR  0.0309274375 2.025976e-05 183 1526.5451 2.001280e-210
#> ETOH 0.0027269529 1.392246e-05 183  195.8671  2.444416e-01
#> GLC  2.2933181627 8.677776e-03 183  264.2749  7.895900e-05
#> LAC  0.0091710408 2.384174e-05 183  384.6633  1.993824e-16
#> k= 3 
#>               rss      var_ref  df     chi2         pval
#> ACE  0.0008564145 3.910221e-06 182 219.0194 3.166424e-02
#> ALA  0.0013111246 6.490101e-06 182 202.0191 1.473762e-01
#> FOR  0.0022969447 1.201390e-05 182 191.1905 3.055006e-01
#> PYR  0.0094076563 2.025976e-05 182 464.3518 1.348948e-26
#> ETOH 0.0026877832 1.392246e-05 182 193.0537 2.733140e-01
#> GLC  2.2180473224 8.677776e-03 182 255.6009 2.623515e-04
#> LAC  0.0044087100 2.384174e-05 182 184.9156 4.258625e-01
#> k= 4 
#>               rss      var_ref  df     chi2         pval
#> ACE  0.0007616405 3.910221e-06 181 194.7819 0.2291696351
#> ALA  0.0012558437 6.490101e-06 181 193.5014 0.2490871435
#> FOR  0.0022975350 1.201390e-05 181 191.2397 0.2866340063
#> PYR  0.0038535634 2.025976e-05 181 190.2078 0.3047174932
#> ETOH 0.0026487379 1.392246e-05 181 190.2492 0.3039794844
#> GLC  2.1327147015 8.677776e-03 181 245.7674 0.0009666041
#> LAC  0.0042670230 2.384174e-05 181 178.9728 0.5286454950
#> k= 5 
#>               rss      var_ref  df chi2      pval
#> ACE  0.0007429421 3.910221e-06 180  190 0.2902766
#> ALA  0.0012331192 6.490101e-06 180  190 0.2902766
#> FOR  0.0022826417 1.201390e-05 180  190 0.2902766
#> PYR  0.0038493543 2.025976e-05 180  190 0.2902766
#> ETOH 0.0026452678 1.392246e-05 180  190 0.2902766
#> GLC  1.6487773464 8.677776e-03 180  190 0.2902766
#> LAC  0.0045299301 2.384174e-05 180  190 0.2902766
#> k= 6 
#>               rss      var_ref  df     chi2       pval
#> ACE  0.0007275983 3.910221e-06 179 186.0760 0.34300006
#> ALA  0.0012226780 6.490101e-06 179 188.3912 0.30048379
#> FOR  0.0022712096 1.201390e-05 179 189.0484 0.28892309
#> PYR  0.0044116953 2.025976e-05 179 217.7565 0.02550377
#> ETOH 0.0026275909 1.392246e-05 179 188.7303 0.29448904
#> GLC  1.7683660757 8.677776e-03 179 203.7810 0.09877058
#> LAC  0.0044048848 2.384174e-05 179 184.7552 0.36838952
```
We can see that $p$-value is above 0.05 for all metabolites only for k=5. We will retain this value as giving the best fit.

The total data fit for k=5 looks like:

``` r
res=resk[[5L]]
nm=colnames(res$mf)[-1L]
matplot(res$tpp, res$isp(res$tpp, nm), t="l")
matpoints(res$tp, res$mf[,nm], pch="o", cex=0.5)
legend("topright", legend=nm, lty=1:5, col=1:6)
```

![plot of chunk unnamed-chunk-24](figure/unnamed-chunk-24-1.png)
and reaction rates:

``` r
matplot(res$tpp, res$vsp(res$tpp), t="l", ylab="Rate [1/s]", xlab="Time")
legend("topright", colnames(bsppar(res$vsp)$qw), lty=1:5, col=1:6, cex=0.75)
```

![plot of chunk unnamed-chunk-25](figure/unnamed-chunk-25-1.png)
The curve oscillations are present for some fluxes but as seen before and based on $\chi^2$ test, in this case, it is not due to overfitting but reflects system behavior through measured data.

In the result files, the rates are plotted with semi-transparent error bands ±2·SD. Let do it here too:

``` r
nmreac=colnames(res$sto)
matplot(res$tpp, res$vsp(res$tpp, nmreac), t="l", ylab="Rate [1/s]", xlab="Time")
tppr=c(res$tpp, rev(res$tpp))
ireac=0
for (reac in nmreac) {
  vp=res$vsp(res$tpp, reac, fsd=2)
  vm=rev(res$vsp(res$tpp, reac, fsd=-2))
  polygon(tppr, c(vp, vm), border=NA, col=do.call(rgb, as.list(c(col2rgb(ireac%%6+1)/255, 0.3))))
  ireac=ireac+1
}
legend("topright", nmreac, lty=1:5, col=1:6, cex=0.75)
```

![plot of chunk unnamed-chunk-26](figure/unnamed-chunk-26-1.png)

Having estimated rates, we can now explore details of metabolic fluxes. Not only we have the total flux $\frac{dM_i}{dt}$ but also all its components coming from involved reactions $S_{ij}v_j$. For PYR, it gives:

``` r
nm="PYR"
jnz=names(which(res$stofull[nm,] != 0)) # non-zero coeffs, i.e. reactions involved in PYR mass balance
fl=t(res$stofull[nm,jnz]*t(res$vsp(res$tpp, jnz))) # each rate v_j is multiplied by S_ij
fl=cbind(res$fsp(res$tpp, nm), fl) # add total flux
matplot(res$tpp, fl, t="l", ylab="Flux [mM/s]", xlab="Time [s]", main="PYR flux components")
legend("topright", legend=c("Total", jnz), lty=1:5, col=1:6)
```

![plot of chunk unnamed-chunk-27](figure/unnamed-chunk-27-1.png)

## Conclusion

In this vignette, we have shown a step-by-step procedure for reaction rate estimation from time-series metabolite measurements and a stoichiometric model. We have seen how package parameters can be chosen, how fit quality can be evaluated and how data can be presented and exploited. Based on this vignette, user should be ready to work with his own data set and stoichiometric model.
