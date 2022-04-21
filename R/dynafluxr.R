#' Function to be called from shell command line
#'
#' @param args Character vector, command line parameters (default
#'   \code{commandArgs(trailingOnly=TRUE)})
#' @details
#'   run \code{cli("-h")} in R or \code{Rscript -e 'dynafluxr::cli()' -h}
#'   in shell to get a help page with available option description
#' @examples
#'   # from shell
#'   # $ Rscript --vanilla -e 'dynafluxr::cli()' -m data_kinetics.tsv -s glycolysis.txt
#'
#'   # from R session
#'   ddir=system.file("data", package="dynafluxr")
#'   meas=file.path(ddir, "data_teusink.tsv")
#'   sto=file.path(ddir, "network_teusink.txt")
#'   res=cli(c("-m", meas, "-s", sto, "--skip", "10"))
#'   tp=res$mf$Time
#'   np=length(tp)
#'   # build tpp, fine time vector for smooth spline plotting
#'   npp=10L # number of plot points in each tp interval
#'   dtp=diff(tp)
#'   tpp=tp[1L]+c(0.,cumsum(rep(diff(tp)/npp, each=npp)))
#'   # plot metabolites
#'   matplot(tpp, res$msp(tpp), type="l")
#'   matpoints(tp, res$mf[,-1], pch="o", cex=0.5)
#'   legend("topright", legend=colnames(bsppar(res$msp)$qw), lty=1:5, col=1:6)
#'   # plot fluxes
#'   dev.new()
#'   matplot(tpp, res$fsp(tpp), type="l")
#'   ref=t(read.table(file.path(ddir, "glyco_teusink.flux.tsv"), sep="\t", header=TRUE, row.names=1, check.names=FALSE))
#'   tf=as.numeric(rownames(ref)) # reference flux time points
#'   nm_flux=colnames(bsppar(res$fsp)$qw)
#'   itf=(tf >= min(tp) & tf <= max(tp))
#'   matpoints(tf[itf], ref[itf, nm_flux, drop=FALSE], pch="o", cex=0.5)
#'   legend("topright", legend=nm_flux, lty=1:5, col=1:6)
#'   # plot residuals
#'   dev.new()
#'   matplot(tpp, res$rsp(tpp), type="l")
#'   legend("topright", legend=colnames(bsppar(res$rsp)$qw), lty=1:5, col=1:6)
#' @importFrom qpdf pdf_combine
#' @importFrom bspline bsppar
#' @importFrom grDevices cairo_pdf
#' @export
cli=function(args=commandArgs(trailingOnly=TRUE)) {
  #stopifnot(requireNamespace("optparse", quite=TRUE))
  library(optparse)
  
  olist=list(
    make_option(c("-m", "--meas"), type="character",
      help="Measurement file in TSV (tab separated value) format,
  one metabolite per column. Time column name should start with 'Time'
  and there should be only one such column.
  Metabolite names should be the same as in stoichiometric model.
  Metabolite measurements whose name is absent in STO file, will be ignored.
  Empty cells are interpreted as NA (non available) and not counted in spline fit.
  If a metabolite has a full column of NA, it will be removed from stoichiometric
  balance and thus from flux least squares. This is different from a metabolite
  which is simply absent from measurement file. In the latter case, a metabolite
  does account in stoichiometric balance but is considered as 0."
    ),
    make_option(c("-s", "--sto"), type="character",
      help="File name with stoichiometric model in plain text format, e.g.:

  vGLK\tGLCi + P -> G6P

  Here 'vGLK' is a reaction name followed by a tab character;
  'GLCi', 'P' and 'G6P' are metabolite names;
  '->' reaction side separator. A symbol '<->' can be used for
  reversible reactions. However, for this application it is
  irrelevant if a reaction is reversible or not;
  '+' is separator of metabolites on the same reaction side.
  Stoichiometric coefficients different from 1 can be used with ' * '
  sign, e.g.:

  vTreha\t2 * G6P + P -> Trh

  here '2' is a such coefficient. The spaces ' ' around '*' are important."
    ),
    make_option(c("-k", "--knot"), type="integer", default=5, help=
      "Internal knot number for B-splines fitting metabolite and flux dynamics
  [default %default]"
    ),
    make_option(c("-n", "--norder"), type="integer", default=4, help=
      "B-splines polynomial order for metabolite kinetics. The flux dynamics will have
  order 'norder-1' [default %default]"
    ),
    make_option(c("-c", "--constr"), type="character",
      help="Constraint file in TSV format, 3 column table: Time, Metabolite, Value.
  Metabolite names should be the same as in stoichiometric model."
    ),
    make_option(c("--mono"), type="character",
      help="Monotonicity file in TSV format, 2 column table: Metabolite, Value 
  (in this order).
  Metabolite names should be the same as in stoichiometric model. 'Value'
  can be 1 (monotonically increasing); 0 (no constraint); -1 (monotonically decreasing).
  For example, substrates that are only consumed could get more
  realistic fit if corresponding 'Value' is set to -1."
    ),
    make_option(c("-o", "--out"), type="character",
      help="Directory (or zip) name to use for result files. By default, measurement name without extension is used."
    ),
    make_option(c("--skip"), type="integer", default=0L,
      help="Number of first time points that should be skipped in metabolite measurements"
    ),
    make_option(c("-z", "--zip"), action="store_true", default=FALSE,
      help="Create zip archive with results."
    ),
    make_option(c("--ils"), action="store_true", default=FALSE,
      help="ils Least Squares formulation"
    ),
    make_option(c("--npp"), type="integer", default=10, help=
      "Number of sub-intervals in each time interval for smooth curve plotting
  [default %default]"
    )
  )
  parser=OptionParser(usage = "Rscript --vanilla -e 'dynafluxr::cli()' -m|--meas MEAS -s|--sto STO [options]
  
  Retrieve flux dynamics from metabolic kinetics. Results are written in a zip file.\n
  Example: Rscript --vanilla -e 'dynafluxr::cli()' -m data_kinetics.tsv -s glycolysis.txt", option_list=olist)
  opt <- try(parse_args(parser, args=args), silent=TRUE)
  if (inherits(opt, "try-error")) {
    if (length(grep("help requested\n$", opt)) == 0L) {
      print_help(parser)
      mes=strsplit(opt, " : ")[[1L]]
      message("Error: ", mes[length(mes)])
    }
    stop(.call=FALSE)
  }
  # mandatory arguments
  if (is.null(opt$meas)) {
    print_help(parser)
    stop("Measurement file must be supplied with '-m' or '--meas' option.", call.=FALSE)
  }
  if (is.null(opt$sto)) {
    print_help(parser)
    stop("Stoichiometric model file must be supplied with '-s' or '--sto' option.", call.=FALSE)
  }
  # get file contents and call fdyn()
  sto=txt2sto(opt$sto) # sto matrix
  mf=read.delim(opt$meas, check.names=FALSE) # metabolic data.frame
  iti=pmatch("Time", colnames(mf))
  if (is.na(iti))
    stop(sprintf("Column name staring with 'Time' is not found or there are many in '%s'", opt$meas))
  #print(c("mf=", colnames(mf)))
  #print(c("sto=", rownames(sto)))
  mf=cbind(Time=mf[,iti], mf[, intersect(colnames(mf)[-iti], rownames(sto))])
  if (ncol(mf) == 1L)
    stop(sprintf("No valid metabolite names in '%s'", opt$meas))
  if (opt$skip > 0L)
    mf=mf[-seq_len(opt$skip),]
  # read equality constraints
  if (!is.null(opt$constr)) {
    dfeq=read.delim(opt$constr, comment.char="#")
    lieq=lapply(colnames(mf)[-1L], function(met) as.matrix(subset(dfeq, Metabolite==met, c(Time, Value))))
    names(lieq)=colnames(mf)[-1L]
  } else {
    lieq=NULL
  }
  # read monotonicity constraints
  if (!is.null(opt$mono)) {
    dfmo=read.delim(opt$mono, comment.char="#")
    mono=structure(double(ncol(mf)-1L), names=colnames(mf)[-1L])
    mono[dfmo[,1L]]=dfmo[,2L]
  } else {
    mono=0
  }
  #print(c("opt=", opt))
  res=fdyn(mf, sto, nsp=opt$norder, nki=opt$knot, lieq=lieq, monotone=mono, ils=opt$ils)
  #res
  # write result files (rd is a temporary dir for results)
  # at the end we'll move all files into a zip archive in the working dir
  #rd=tempfile(pattern="fdyn")
  if (is.null(opt$out)) {
    rd=tools::file_path_sans_ext(opt$meas)
  } else {
    rd=opt$out
  }
  if (!dir.exists(rd))
    dir.create(rd)
  #bnm=tools::file_path_sans_ext(opt$meas)
  pm=bspline::bsppar(res$msp) # metab params
  pf=bspline::bsppar(res$fsp) # flux params
  pr=bspline::bsppar(res$rsp) # resid params
  pi=bspline::bsppar(res$isp) # resid params
  tp=mf[,1L]
  tpp=tp[1L]+c(0.,cumsum(rep(diff(tp)/opt$npp, each=opt$npp)))
  ratp=range(tp, na.rm=TRUE)
  # pdf with metabs
  pdf(file.path(rd, "met.pdf"))
  mc=res$msp(tpp) # metabolite smooth curves
  plot(1, xlim=ratp, main="Measured concentrations fitted by B-splines", ylim=range(mc, mf[,-1L], na.rm=TRUE), xlab="Time", ylab="Concentration", t="n")
  matlines(tpp, mc)
  matpoints(tp, res$mf[,-1], pch="o", cex=0.5)
  legend("topright", legend=colnames(pm$qw), lty=1:5, col=1:6)
  for (m in colnames(pm$qw)) {
    if (all(is.na(mf[,m])))
      next
    plot(1, main=m, xlab="Time", ylab="Concentration", xlim=ratp, ylim=range(mf[,m], mc[,m], na.rm=TRUE), type="n")
    points(tp, mf[,m], pch="o", cex=0.5)
    lines(tpp, mc[,m])
  }
  dev.off()
  # pdf with fluxes
  pdf(file.path(rd, "flux.pdf"))
  fc=res$fsp(tpp) # flux smooth curves
  matplot(tpp, fc, xlab="Time", ylab="Flux", type="l")
  legend("topright", legend=colnames(pf$qw), lty=1:5, col=1:6)
  for (f in colnames(pf$qw)) {
    plot(1, main=f, xlab="Time", ylab="Flux", xlim=ratp, ylim=range(fc[,f], na.rm=TRUE), type="n")
    lines(tpp, fc[,f])
  }
  dev.off()
  # pdf with restored (integrated) metabs
  grDevices::cairo_pdf(file.path(rd, "imet%03d.pdf"))
  ic=res$isp(tpp) # imet smooth curves
  plot(1, main="Estimated concentrations", xlim=ratp, ylim=range(ic, mf[,-1L], na.rm=TRUE), xlab="Time", ylab="\u222bS\u00b7f dt", t="n")
  matlines(tpp, ic)
  tmp=matrix(NA_real_, nrow(mf), ncol(ic)) # inject here mf values
  colnames(tmp)=colnames(ic)
  cnm=intersect(colnames(mf)[-1L], colnames(ic)) # NA columns in mf are absent in ic
  tmp[,cnm]=as.matrix(mf[,cnm])
  matpoints(tp, tmp, pch="o", cex=0.5)
  legend("topright", legend=colnames(pr$qw), lty=1:5, col=1:6)
  for (m in colnames(pi$qw)) {
    plot(1, main=m, xlab="Time", ylab="\u222bS\u00b7f dt", xlim=ratp, ylim=range(ic[,m], tmp[,m], na.rm=TRUE), type="n")
    lines(tpp, ic[, m])
    points(tp, tmp[, m], pch="o", cex=0.55)
  }
  dev.off()
  li=list.files(rd, "imet[0-9]+\\.pdf", full.names=TRUE)
  qpdf::pdf_combine(input=li, output=file.path(rd, "imet.pdf"))
  unlink(li)
  # pdf with residuals
  pdf(file.path(rd, "resid.pdf"))
  rc=res$rsp(tpp) # residual smooth curves
  matplot(tpp, rc, xlab="Time", ylab="dm/dt - S\u00b7f", type="l")
  legend("topright", legend=colnames(pr$qw), lty=1:5, col=1:6)
  for (r in colnames(pr$qw)) {
    plot(1, main=r, xlab="Time", ylab="dm/dt - S\u00b7f", xlim=ratp, ylim=range(rc[,r], na.rm=TRUE), type="n")
    lines(tpp, rc[,r])
  }
  dev.off()
  # tsv
  write.table(cbind(Time=tpp, mc), file=file.path(rd, "met.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(cbind(Time=tpp, fc), file=file.path(rd, "flux.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  # .RData
  save(res, file=file.path(rd, "env.RData"))
  # zip files
  if (opt$zip) {
    zip(paste0(rd, ".zip"), rd, extra="-j")
    unlink(rd, recursive=TRUE)
  }
  invisible(res)
}

#' Retrieve flux dynamics from metabolic kinetics
#' 
#' @param mf Data-frame or matrix, metabolite kinetic measurements.
#'   Columns must be named with metabolite names and 'Time'.
#' @param sto Stoichiometric matrix, \code{sto[i,j]} means reaction 'j' produces
#'   metabolite 'i' with the flux 'sto[i,j]'. If \code{sto[i,j] < 0},
#'   the metabolite 'i' is consumed. Columns must be named with flux names.
#'   Rows must be names with metabolite names.
#' @param nsp Integer, polynomial order of B-spline to use for metabolites
#' @param nki Integer, number of internal knots for B-splines
#' @param lieq List, equality constraints on metabolites
#' @param monotone Numeric scalar or vector, 1=metabolites are
#'   monotonically increasing;
#'   -1=monotonically decreasing; 0=no constraint. If vector, each value
#'   constraints (or not) a corresponding data column in mf ('Time'
#'   column is excluded
#'   from counting)
#' @details
#'   Each item in \code{lieq} corresponds to a metabolite and is a
#'   2 column matrix (Time, Value). Each
#'   row of this matrix indicates what 'Value' must take corresponding
#'   metabolite at what 'Time'. Typically, it can be used to impose
#'   starting values at Time=0 for some metabolites.\cr
#'   All metabolite fits are constraint to have values >= 0.
#' @return List with following components:
#' \describe{
#'   \item{mf:}{ metabolite data frame used for fitting}
#'   \item{sto:}{ stoichiometric matrix used for fitting}
#'   \item{msp:}{ metabolite spline function}
#'   \item{fsp:}{ flux spline function}
#'   \item{dsp:}{ metabolite first derivative spline function}
#'   \item{rsp:}{ residual \code{dm/dt - sto\%*\%flux} spline function}
#' }
#' @importFrom bspline fitsmbsp dbsp bsppar par2bsp
#' @export
fdyn=function(mf, sto, nsp=4L, nki=5L, lieq=NULL, monotone=0, ils=FALSE) {

  tp=mf$Time
  dtp=diff(tp)
  np=length(tp)
  # fit measured metabs
  msp=bspline::fitsmbsp(tp, mf[, -1L, drop=FALSE], n=nsp, nki=nki, monotone=monotone, positive=1, lieq=lieq, control=list(monotone=TRUE, errx=dtp[1]/10., trace=1))
  # remove metab's NA
  ina=names(which(apply(bspline::bsppar(msp)$qw, 2, function(vc) anyNA(vc))))
  if (length(ina)) {
    sto=sto[-match(ina, rownames(sto)),,drop=FALSE]
    mf=mf[,-match(ina, colnames(mf)),drop=FALSE]
    e=environment(msp)
    e$qw=e$qw[,-match(ina, colnames(e$qw)),drop=FALSE]
  }
  parm=bspline::bsppar(msp)
  nmet=nrow(sto)
  nflux=ncol(sto)
  nwm=nrow(parm$qw)
  # first derivatives
  dsp=bspline::dbsp(msp)
  pard=bspline::bsppar(dsp)
  qwd=pard$qw
  nwf=nrow(pard$qw)
  # complete derivatives by 0 for non measured metabs
  qwd0=matrix(0., nrow(qwd), nmet)
  colnames(qwd0)=rownames(sto)
  qwd0[,colnames(qwd)]=qwd

  # calculates weights for fluxes
  if (ils) {
    # build integral LS
    # qwm0: metab coeffs including 0s
    qwm0=matrix(0., nrow=nwm, ncol=nmet)
    colnames(qwm0)=rownames(sto)
    qwm0[,colnames(parm$qw)]=parm$qw
    # unknowns: qwf and mstart
    mstart=setNames(double(nmet), colnames(qwm0))
    mstart[colnames(parm$qw)]=parm$qw[1,] # starting values where known
    # integration matrix (cumsum(...))
    ima=diag(nrow=nwm-1)
    ima=(row(ima)>=col(ima))+0.
    ima=rbind(0., arrApply::arrApply(ima, 2, "multv", v=diff(pard$xk, lag=nsp)/nsp))
    # f2m: ima%*%qwf%*%t(sto)+rep(1,nkmet)%o%mstart
    jqw=aperm(ima%o%t(sto), c(1L,4L,2L,3L))
    dim(jqw)=c(prod(dim(jqw)[1L:2L]), prod(dim(jqw)[3L:4L]))
    jadd=rep(1., nrow(ima))%o%diag(nrow=nmet)
    dim(jadd)=c(nrow(ima)*nmet, nmet)
    jtot=cbind(jadd, jqw)
    # inequalities: qwm >= 0; dm/dt monotone for some metabs
    if (any(monotone != 0)) {
      if (length(monotone) == 1L)
        monotone=setNames(rep(monotone, nmet), rownames(sto))
      i=names(which(monotone > 0))
      umono=aperm(diag(nrow=nwf)%o%sto[i,,drop=FALSE], c(1L,3L,2L,4L))
      dim(umono)=c(nwf*length(i), nwf*nflux)
      i=names(which(monotone < 0))
      tmp=aperm(diag(nrow=nwf)%o%(-sto[i,,drop=FALSE]), c(1L,3L,2L,4L))
      dim(tmp)=c(nwf*length(i), nwf*nflux)
      umono=cbind(matrix(0., nrow=nrow(umono)+nrow(tmp), ncol=nmet), rbind(umono, tmp))
      cmono=double(nrow(umono))
    } else {
      umono=NULL
      cmono=NULL
    }
    # positivity matrix u is the same as jtot
    u=rbind(jtot, umono)
    co=double(nrow(u))
    # solve iLS
    st=system.time({p=nlsic::lsi(jtot, c(qwm0), u=u, co=co)})
    icnst=seq_len(nmet)
    mst=p[icnst]
    qwf=p[-icnst]
    dim(qwf)=c(nwf, nflux)
    colnames(qwf)=colnames(sto)
  } else {
    # differential LS
    # generalized inverse
    st=system.time({
    s=svd(sto)
    d=s$d
    d[d <= d[1L]*1.e-10]=0.
    d[d != 0.]=1./d[d != 0.]
    stoinv=s$v%*%(d*t(s$u))
    dimnames(stoinv)=rev(dimnames(sto))
    qwf=qwd0%*%t(stoinv)
    })
  }
  #print(st)
  # flux splines
  fsp=bspline::par2bsp(nsp-1L, qwf, pard$xk)
  # metab restored from S*f
  qwde=qwf%*%t(sto) # dm/dt estimated from fluxes
  if (ils) {
    const=mst
  } else {
    const=setNames(double(nmet), colnames(qwde))
    const[colnames(parm$qw)]=parm$qw[1,] # starting values where known
  }
  isp=ibsp(bspline::par2bsp(nsp-1L, qwde, pard$xk), const=const)
  # residuals dm/dt-sto*f
  rsp=bspline::par2bsp(nsp-1L, qwd0-qwde, pard$xk)
  list(mf=mf, sto=sto, invsto=if (ils) NULL else stoinv, msp=msp, fsp=fsp, dsp=dsp, isp=isp, rsp=rsp)
}
