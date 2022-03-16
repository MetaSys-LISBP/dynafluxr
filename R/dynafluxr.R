#' Function to be called from shell command line
#'
#' @examples
#'   ddir=system.file("data", package="dynafluxr")
#'   meas=file.path(ddir, "data_teusink.tsv")
#'   sto=file.path(ddir, "network_teusink.txt")
#'   res=cli(c("-m", meas, "-s", sto, "--skip", "10"))
#'   tp=res$mf$Time
#'   np=length(tp)
#'   npp=10L # number of plot points in each tp interval
#'   dtp=diff(tp)
#'   tpp=tp[1L]+c(0.,cumsum(rep(diff(tp)/npp, each=npp)))
#'   # plot metabolites
#'   matplot(tpp, res$msp(tpp), type="l")
#'   matpoints(tp, res$mf[,-1], pch=".")
#'   legend("topright", legend=colnames(bsppar(res$msp)$qw), lty=1:5, col=1:6)
#'   # plot fluxes
#'   dev.new()
#'   matplot(tpp, res$fsp(tpp), type="l")
#'   ref=t(read.table(file.path(ddir, "glyco_teusink.flux.tsv"), sep="\t", header=TRUE, row.names=1, check.names=FALSE))
#'   tf=as.numeric(rownames(ref)) # reference flux time points
#'   nm_flux=colnames(bsppar(res$fsp)$qw)
#'   itf=(tf >= min(tp) & tf <= max(tp))
#'   matpoints(tf[itf], ref[itf, nm_flux, drop=FALSE], pch=".")
#'   legend("topright", legend=nm_flux, lty=1:5, col=1:6)
#'   # plot residuals
#'   dev.new()
#'   matplot(tpp, res$rsp(tpp), type="l")
#'   legend("topright", legend=colnames(bsppar(res$rsp)$qw), lty=1:5, col=1:6)
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
  Metabolite measurements whose name is absent in STO file, will be ignored."
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
     make_option(c("--skip"), type="integer", default=0L,
      help="Number of first time points that should be skipped in metabolite measurements"
    )
  )
  parser=OptionParser(usage = "Rscript --vanilla -e 'dynafluxr::cli()' -m|--meas MEAS -s|--sto STO [options]
  
  Retrieve flux dynamics from metabolic kinetics\n
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
  } else {
    lieq=NULL
  }
  #print(c("opt=", opt))
  res=fdyn(mf, sto, nsp=opt$norder, nki=opt$knot, lieq)
  res
  # write result files
}

#' Retrieve flux dynamics from metabolic kinetics
#' 
#' @return List with names components: \itemize{
#'   \item{\code{mf}}{: metabolite data frame used for fitting}
#'   \item{\code{sto}}{: stoichiometric matrix used for fitting}
#'   \item{\code{msp}}{: metabolite spline function}
#'   \item{\code{fsp}}{: flux spline function}
#'   \item{\code{dsp}}{: metabolite first derivative spline function}
#'   \item{\code{rsp}}{: residual \code{dm/dt - sto%*%flux} spline function}
#' }
#' @export

fdyn=function(mf, sto, nsp=4L, nki=5L, lieq=NULL) {
  nmet=nrow(sto)
  nflux=ncol(sto)

  tp=mf$Time
  dtp=diff(tp)
  np=length(tp)
  # generalized inverse
  s=svd(sto)
  d=s$d
  d[d <= d[1]*1.e-10]=0.
  d[d != 0.]=1./d[d != 0.]
  stoinv=s$v%*%(d*t(s$u))
  dimnames(stoinv)=rev(dimnames(sto))

  msp=fitsmbsp(tp, mf[, -1L, drop=FALSE], n=nsp, control=list(monotone=TRUE, errx=dtp[1]/10.), lieq=lieq)
  # first derivatives
  dsp=dbsp(msp)
  pard=bsppar(dsp)
  qwd=pard$qw
  # complete mqw2 by 0 for non measured metabs
  mfqwd=matrix(0., nrow(qwd), nmet) # 'f' for full
  colnames(mfqwd)=rownames(sto)
  mfqwd[,colnames(qwd)]=qwd

  # calculates weights for fluxes
  qwf=mfqwd%*%t(stoinv)
  # flux splines
  fsp=par2bsp(nsp-1L, qwf, pard$xk)
  # residuals dm/dt-sto*f
  rsp=par2bsp(nsp-1L, mfqwd-qwf%*%t(sto), pard$xk)
  list(mf=mf, sto=sto, msp=msp, fsp=fsp, dsp=dsp, rsp=rsp)
}
