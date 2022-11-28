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
#'   ddir=system.file("dataglyco", package="dynafluxr")
#'   meas=file.path(ddir, "data_teusink.tsv")
#'   sto=file.path(ddir, "network_teusink.txt")
#'   res=cli(c("-m", meas, "-s", sto, "--skip", "10", "-o", ""))
#'   tp=res$tp
#'   np=length(tp)
#'   tpp=res$tpp
#'   # plot species
#'   matplot(tpp, res$msp(tpp), type="l")
#'   matpoints(tp, res$mf[,-1], pch="o", cex=0.5)
#'   legend("topright", legend=colnames(bsppar(res$msp)$qw), lty=1:5, col=1:6, cex=0.75)
#'   # plot rates
#'   dev.new()
#'   matplot(tpp, res$vsp(tpp), type="l")
#'   ref=t(read.delim(file.path(ddir, "glyco_teusink.flux.tsv"), row.names=1, check.names=FALSE))
#'   tf=as.numeric(rownames(ref)) # reference rate time points
#'   nm_rate=colnames(bsppar(res$vsp)$qw)
#'   itf=(tf >= min(tp) & tf <= max(tp))
#'   matpoints(tf[itf], ref[itf, nm_rate, drop=FALSE], pch="o", cex=0.5)
#'   legend("topright", legend=nm_rate, lty=1:5, col=1:6, cex=0.75)
#'   # plot residuals
#'   dev.new()
#'   matplot(tpp, res$risp(tpp), type="l")
#'   legend("topright", legend=colnames(bsppar(res$rsp)$qw), lty=1:5, col=1:6, cex=0.75)
#' @importFrom qpdf pdf_combine
#' @importFrom bspline bsppar
#' @importFrom grDevices cairo_pdf rgb col2rgb pdf dev.off
#' @importFrom optparse make_option OptionParser parse_args print_help
#' @importFrom utils read.delim write.table zip
#' @importFrom stats setNames
#' @importFrom graphics legend matlines matpoints polygon points lines 
#' @export
cli=function(args=commandArgs(trailingOnly=TRUE)) {
  Specie <- Time <- Value <- NULL # to keep 'R CMD check' calm about subset() arguments
  olist=list(
    make_option(c("-m", "--meas"), type="character",
      help="Measurement file in TSV (tab separated value) format,
  one specie per column. Time column name should start with 'Time'
  and there should be only one such column.
  Specie names should be the same as in stoichiometric model.
  Specie measurements whose name is absent in STO file, will be ignored.
  Empty cells are interpreted as NA (non available) and not counted in spline fit.
  If a specie has a full column of NA, it will be removed from stoichiometric
  balance and thus from rate least squares. This is different from a specie
  which is simply absent from measurement file. In the latter case, a specie
  does account in stoichiometric balance but is considered as 0."
    ),
    make_option(c("-s", "--sto"), type="character",
      help="File name with stoichiometric model in plain text format, e.g.:

  vGLK\tGLCi + P -> G6P

  Here 'vGLK' is a reaction name followed by a tab character;
  'GLCi', 'P' and 'G6P' are specie names;
  '->' reaction side separator. A symbol '<->' can be used for
  reversible reactions. However, for this application it is
  irrelevant if a reaction is reversible or not;
  '+' is separator of species on the same reaction side.
  Stoichiometric coefficients different from 1 can be used with ' * '
  sign, e.g.:

  vTreha\t2 * G6P + P -> Trh

  here '2' is a such coefficient. The spaces ' ' around '*' are important."
    ),
    make_option(c("-k", "--knot"), type="integer", default=5, help=
      "Internal knot number for B-splines fitting specie and rate dynamics
  [default %default]"
    ),
    make_option(c("-n", "--norder"), type="integer", default=4, help=
      "B-splines polynomial order for specie kinetics. The rate dynamics will have
  order 'norder-1' [default %default]"
    ),
    make_option(c("--lna"), type="character", default="",
      help="List of coma separated NA species, e.g. '--lna=FBP,F6P'. If a specie is present both in MEAS file and in LNA list, it is overwritten with NAs"
    ),
    make_option(c("-c", "--constr"), type="character",
      help="Constraint file in TSV format, 3 column table: Time, Specie, Value.
  Specie names should be the same as in stoichiometric model."
    ),
    make_option(c("-a", "--atom"), type="character",
      help="Atom length file in TSV format, 2 column table: Specie, Atom_length (in this order).
  Specie names should be the same as in stoichiometric model. Atom_length
  column must contain integer non negative values."
    ),
    make_option(c("--mono"), type="character",
      help="Monotonicity file in TSV format, 2 column table: Specie, Value 
  (in this order).
  Specie names should be the same as in stoichiometric model. 'Value'
  can be 1 (monotonically increasing); 0 (no constraint); -1 (monotonically decreasing).
  For example, substrates that are only consumed could get more
  realistic fit if corresponding 'Value' is set to -1. Cf. also options --increasing and --decreasing."
    ),
    make_option(c("--increasing"), type="character", default="",
      help="List of coma separated species that are supposed to be monotonously increasing, e.g. '--increasing=LAC,ETOH'. If a specie is present both in file MONO and in this option,
  this option takes the precedence."
    ),
    make_option(c("--decreasing"), type="character", default="",
      help="List of coma separated species that are supposed to be monotonously decreasing, e.g. '--decreasing=GLC'. If a specie is present both in file MONO and in this option,
  this option takes the precedence."
    ),
    make_option(c("-o", "--out"), type="character",
      help="Directory (or zip) name to use for result files. By default, measurement name without extension is used. If empty, no results are written to disk (can be useful for programmatic use)."
    ),
    make_option(c("--skip"), type="integer", default=0L,
      help="Number of first time points that should be skipped in specie measurements"
    ),
    make_option(c("-z", "--zip"), action="store_true", default=FALSE,
      help="Create zip archive with results (default: FALSE)."
    ),
    make_option(c("--dls"), action="store_true", default=FALSE,
      help="use Differential Least Squares formulation (default: FALSE)"
    ),
    make_option(c("--npp"), type="integer", default=10, help=
      "Number of sub-intervals in each time interval for smooth curve plotting
  [default %default]"
    ),
    make_option(c("--fsd"), type="double", default=2., help=
      "SD factor for plotting gray band \u00b1fsd*SD around spline curves. Use '--fsd=0' to cancel these bands. [default %default]"
    )
  )
  parser=OptionParser(usage = "Rscript --vanilla -e 'dynafluxr::cli()' -m|--meas MEAS -s|--sto STO [options]
  
  Retrieve rate dynamics from metabolic kinetics. Results are written in a directory or zip file.\n
  Example: Rscript --vanilla -e 'dynafluxr::cli()' -m data_kinetics.tsv -s glycolysis.txt", option_list=olist)
  opt <- try(parse_args(parser, args=args), silent=TRUE)
  if (inherits(opt, "try-error")) {
    if (length(grep("help requested\n$", opt)) == 0L) {
      print_help(parser)
      mes=strsplit(opt, " : ")[[1L]]
      message("Error: ", mes[length(mes)])
    }
    #stop(.call=FALSE)
    if (!base::interactive()) {
      q("no", status=1)
    } else {
      #options(show.error.messages = FALSE)
      #stop("\r", call.=FALSE)
      return(invisible(NULL))
    }
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
  # put Time in first column
  mf=cbind(Time=mf[,iti], mf[, intersect(colnames(mf)[-iti], rownames(sto))])
  if (ncol(mf) == 1L)
    stop(sprintf("No valid specie names in '%s'", opt$meas))
  if (opt$skip > 0L)
    mf=mf[-seq_len(opt$skip),]
  # read equality constraints
  if (!is.null(opt$constr)) {
    dfeq=read.delim(opt$constr, comment.char="#")
    lieq=lapply(colnames(mf)[-1L], function(met) as.matrix(subset(dfeq, Specie==met, c(Time, Value))))
    names(lieq)=colnames(mf)[-1L]
  } else {
    lieq=NULL
  }
  # read monotonicity constraints
  mono=setNames(double(nrow(sto)), rownames(sto))
  if (!is.null(opt$mono)) {
    dfmo=read.delim(opt$mono, comment.char="#")
    mono[dfmo[,1L]]=dfmo[,2L]
  }
  mono[strsplit(opt$increasing, ",")[[1L]]]=1
  mono[strsplit(opt$decreasing, ",")[[1L]]]=-1
  # read atom lengths
  if (!is.null(opt$atom)) {
    datom=read.delim(opt$atom, comment.char="#")
    atomlen=setNames(double(nrow(sto)), rownames(sto))
    atomlen[datom[,1L]]=datom[,2L]
  } else {
    atomlen=NULL
  }
  # prepare NA list form --lna
  lna=strsplit(opt$lna, ",")[[1L]]
  if (length(lna) > 0L) {
    mf[,lna]=rep(NA, nrow(mf))
  }
  #print(c("opt=", opt))
  res=fdyn(mf, sto, nsp=opt$norder, nki=opt$knot, lieq=lieq, monotone=mono, dls=opt$dls, atomlen=atomlen, npp=opt$npp)
  #res
  # write result files (rd is a temporary dir for results)
  # at the end we'll move all files into a zip archive in the working dir
  #rd=tempfile(pattern="fdyn")
  if (is.null(opt$out)) {
    rd=tools::file_path_sans_ext(opt$meas)
  } else {
    rd=opt$out
  }
  if (nchar(rd)) {
    if (!dir.exists(rd))
      dir.create(rd)
    tp=res$tp
    tpp=res$tpp
    ratp=range(tp, na.rm=TRUE)
    tppr=c(tpp, rev(tpp)) # for SD band plotting
    bcol=do.call(rgb, as.list(c(col2rgb(1)/255, 0.3))) # band color
    # define plot function
    plotsp=function(fname, sp, main, ylab, data=NULL, sto=NULL) {
      # plot splines with sd-band and data
      open_here=!is.null(fname)
      if (open_here)
        pdf(fname)
      p=bspline::bsppar(sp)
      mc=sp(tpp) # specie smooth curves
      ylim=range(mc, data, na.rm=TRUE)
      if (!is.null(p$sdqw)) {
        msdp=sp(tpp, fsd=opt$fsd)
        msdm=sp(tpp, fsd=-opt$fsd)
        ylim=range(ylim, msdp, msdm)
      }
      plot(1, xlim=ratp, main=main, ylim=ylim, xlab="Time", ylab=ylab, t="n")
      matlines(tpp, mc)
      if (!is.null(data))
        matpoints(tp, data, pch="o", cex=0.5)
      legend("topright", legend=colnames(p$qw), lty=1:5, col=1:6, bg=rgb(1,1,1,0.3), cex=0.75)
      # multi-color sd-bands
      if (!is.null(p$sdqw)) {
        for (im in seq_along(colnames(p$qw))) {
          m=colnames(p$qw)[im]
          #cat("m=", m, "\n")
          polygon(tppr, c(msdp[,m], rev(msdm[,m])), border=NA, col=do.call(rgb, as.list(c(col2rgb((im-1)%%6+1)/255, 0.3))))
        }
      }
      # individual data + sd-bands
      vct=t(res$vsp(tpp)) # transposed rate smooth curves
      for (m in colnames(p$qw)) {
        if (!(is.null(data) || !(m %in% colnames(data)) || all(is.na(data[,m])))) {
          d=data[,m]
        } else {
          d=NULL
        }
        if (!is.null(sto)) {
          # prepare individual fluxes
          nm_re=names(which(sto[m,] != 0))
          fl=t(sto[m,nm_re]*vct[nm_re,,drop=FALSE])
        } else {
          fl=NULL
        }
        ylim=range(d, mc[,m], fl, na.rm=TRUE)
        if (!is.null(p$sdqw))
          ylim=range(ylim, msdp[,m], msdm[,m])
        plot(1, main=m, xlab="Time", ylab=ylab, xlim=ratp, ylim=ylim, type="n")
        if (!is.null(d))
          points(tp, d, pch="o", cex=0.5)
        lines(tpp, mc[,m])
        if (!is.null(p$sdqw))
          polygon(tppr, c(msdp[,m], rev(msdm[,m])), border=NA, col=bcol)
        if (!is.null(sto)) {
          matlines(tpp, fl, lty=c(2:5,1), col=c(2:6,1))
          legend("topright", legend=c("Total", nm_re), lty=1:5, col=1:6)
        }
      }
      if (open_here)
        dev.off()
      return(invisible(mc))
    }
    # pdf with species
    mc=plotsp(file.path(rd, "specie.pdf"), res$msp, "Measured concentrations fitted by B-splines", "Concentration", res$mf[,-1L])
    # pdf with atom balance
    if (length(atomlen) > 0L) {
      pdf(file.path(rd, "atom.pdf"))
      datom=colSums(t(mf[,-1L])*atomlen[colnames(mf)[-1L]], na.rm=TRUE)
      ac=res$asp(tpp)
      iac=res$iasp(tpp)
      plot(1, xlim=ratp, main="Atom balance evolution", ylim=range(ac, iac, datom, na.rm=TRUE), xlab="Time", ylab="Total atom number", t="n")
      matlines(tpp, cbind(ac, iac))
      points(tp, datom, pch="o", cex=0.5)
      legend("topright", legend=c("fitted species", "integrated species"), bg=rgb(1,1,1,0.3), lty=1:2, col=1:2, cex=0.75)
      dev.off()
      atomline=" - `atom.pdf`: atom balance plots;"
    } else {
      atomline=NULL
    }
    # pdf with rates
    #browser()
    vc=plotsp(file.path(rd, "rate.pdf"), res$vsp, "Reaction rates", "Rate 1/[Time]", NULL)
    # pdf with fluxes (S*v)
    fc=plotsp(file.path(rd, "flux.pdf"), res$fsp, "Total fluxes (S\u00b7v)", "Flux [Conctr]/[Time]", NULL, res$stofull)
    # pdf with restored (integrated) metabs
    grDevices::cairo_pdf(file.path(rd, "imet%03d.pdf"))
    inames=colnames(bspline::bsppar(res$isp)$qw)
    tmp=matrix(NA_real_, nrow(mf), length(inames)) # inject here mf values
    colnames(tmp)=inames
    cnm=intersect(inames, colnames(mf)[-1L]) # NA columns in mf are absent in mc
    tmp[,cnm]=as.matrix(mf[,cnm])
    ic=plotsp(NULL, res$isp, "Estimated concentrations", "\u222bS\u00b7v dt", tmp)
    dev.off()
    li=list.files(rd, "imet[0-9]+\\.pdf", full.names=TRUE)
    qpdf::pdf_combine(input=li, output=file.path(rd, "ispecie.pdf"))
    unlink(li)
    # pdf with residuals
    rc=plotsp(file.path(rd, "resid.pdf"), res$rsp, "Residuals", "dm/dt - S\u00b7v", NULL)
    # tsv
    write.table(cbind(Time=tpp, mc), file=file.path(rd, "specie.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
    write.table(cbind(Time=tpp, ic), file=file.path(rd, "ispecie.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
    write.table(cbind(Time=tpp, vc), file=file.path(rd, "rate.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
    write.table(cbind(Time=tpp, fc), file=file.path(rd, "flux.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
    # .RData
    save(res, file=file.path(rd, "env.RData"))
    # Readme.md
    cat(file=file.path(rd, "Readme.md"), sep="\n",
      "# Retrieving reaction rate dynamics from specie kinetics (dynafluxr results)", "",
      paste0("This is the result files produced by dynafluxr R package on ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z (%Z).")), "",
      "The command to reproduce these results is:", "",
      paste0("`Rscript --vanilla -e 'dynafluxr::cli()' ", paste0(shQuote(args), collapse=" "), "`"),
      "", "##File contents", "",
      " - `specie.pdf`: concentration plots (fitted by B-spline);",
      " - `ispecie.pdf`: estimated concentration plots vs Time (by integration of *S\u00b7v*);",
      atomline,
      " - `rate.pdf`: estimated rate plots (by solving least squares);",
      " - `flux.pdf`: estimated total fluxes (S*v) plots (by solving least squares);",
      " - `resid.pdf`: residuals *dm/dt - S\u00b7v* plots;",
      " - `specie.tsv`: concentration table;",
      " - `ispeci.tsv`: estimated concentration table;",
      " - `rate.tsv`: rate table;",
      " - `flux.tsv`: flux table;",
      " - `env.RData`: stored R list `res` such as returned by `dynafluxr::fdyn()`. It can be read in R session with `e=new.env(); load('env.RData', envir=e)` and then accessed as e.g. `ls(e$res)`;",
      " - `Readme.md`: this file;"
    )
    # zip files
    if (opt$zip) {
      zip(paste0(rd, ".zip"), rd, extras="-j")
      unlink(rd, recursive=TRUE)
    }
  }
  invisible(res)
}

#' Retrieve flux dynamics from metabolic kinetics
#' 
#' @param mf Data-frame or matrix, specie kinetic measurements.
#'   Columns must be named with specie names and 'Time'.
#' @param stofull Full stoichiometric matrix, \code{stofull[i,j]} means reaction 'j' produces
#'   specie 'i' with coeficient 'stofull[i,j]'. If \code{stofull[i,j] < 0},
#'   the specie 'i' is consumed. Columns must be named with reaction names.
#'   Rows must be names with the species. "Full" in the name means that matrix includes
#'   even NA species.
#' @param nsp Integer, polynomial order of B-spline to use for species
#'   (default 4)
#' @param nki Integer, number of internal knots for B-splines
#'   (default 5)
#' @param lieq List, equality constraints on species
#'   (default NULL, i.e. no equality constraint)
#' @param monotone Numeric scalar or vector, 1=species are
#'   monotonically increasing;
#'   -1=monotonically decreasing; 0=no constraint. If vector, each value
#'   constraints (or not) a corresponding data column in mf ('Time'
#'   column is excluded
#'   from counting)
#'   (default 0, i.e. no monotonicity constraint)
#' @param dls Logical scalar, if TRUE, indicates that differential least squares
#'   should be resolved instead of integral least squares.
#'   (default FALSE, i.e. ILS will be used)
#' @param atomlen Numerical named vector, indicates what is label length
#'   of a given specie used a vector item name. If provided, results
#'   will contain \code{lsp} and \code{ilsp} fields which
#'   are a B-spline function representing atom balance over msp and isp splines.
#'   (default NULL, i.e. no atom balance will be provided)
#' @param npp Integer scalar, indicates a number of subintervals to which
#'   each time interval will be subdivided to produce smoother plots.
#'   (default 10)
#' @details
#'   Each item in \code{lieq} corresponds to a specie and is a
#'   2 column matrix (Time, Value). Each
#'   row of this matrix indicates what 'Value' must take corresponding
#'   specie at what 'Time'. Typically, it can be used to impose
#'   starting values at Time=0 for some species.\cr
#'   All specie fits are constraint to have values >= 0.
#' @return List with following components:
#' \describe{
#'   \item{mf:}{ specie data frame used for fitting}
#'   \item{tp:}{ vector of time points for used measurements}
#'   \item{tpp:}{ vector of time points for plot (fine time resolution)}
#'   \item{sto:}{ stoichiometric matrix used for fitting}
#'   \item{invsto:}{ pseudo-inverse of stoichiometric matrix used in DSL (NULL for ISL)}
#'   \item{msp:}{ measured specie spline function}
#'   \item{vsp:}{ estimated rates spline function}
#'   \item{fsp:}{ estimated total flux (S*v) spline function}
#'   \item{dsp:}{ first derivative of measured spline function}
#'   \item{isp:}{ integrated specie spline function}
#'   \item{asp:}{ atom balance over msp spline function}
#'   \item{iasp:}{ atom balance over isp spline function}
#'   \item{vsp:}{ flux spline function}
#'   \item{dsp:}{ measured specie first derivative spline function}
#'   \item{rsp:}{ residual \code{dM/dt - S\%*\%v} spline function}
#'   \item{risp:}{ integral residual \code{M - \\u222bS\%*\%v dt} spline function}
#'   \item{sdrate:}{ matrix of SD values for flux B-spline coefficients, of size (\code{ncoef x nrate})}
#'   \item{chi2tab:}{ data-frame with chi2-test results}
#' }
#' @importFrom bspline fitsmbsp dbsp bsppar par2bsp ibsp
#' @importFrom slam simple_triplet_zero_matrix matprod_simple_triplet_matrix
#' @importFrom arrApply arrApply
#' @importFrom stats var pchisq
#' @importFrom nlsic lsi
#' @export
fdyn=function(mf, stofull, nsp=4L, nki=5L, lieq=NULL, monotone=0, dls=FALSE, atomlen=NULL, npp=10L) {

  tp=mf$Time
  dtp=diff(tp)
  np=length(tp)
  tpp=tp[1L]+c(0.,cumsum(rep(dtp/npp, each=npp)))
  # prepare mono
  if (length(monotone) > 1L) {
    mono=monotone[colnames(mf)[-1L]]
  } else {
    mono=monotone
  }
  # estimate var_ref with biggest d2 in variance (for chi2 test)
  nki_test=seq(max(0, nki-5), nki+5)
  var_test=sapply(nki_test, function(k) {s=bspline::smbsp(tp, mf[, -1L, drop=FALSE], n=nsp, nki=k, monotone=mono, positive=1, lieq=lieq, estSD=FALSE); colMeans((s(tp)-mf[, -1L, drop=FALSE])**2, na.rm=TRUE)})
  i_ref=which.max(base::diff(base::colSums(var_test, na.rm=TRUE), difference=2L))+1L
  var_ref=na.omit(var_test[,i_ref])
  #print(list(nki_test, var_test, i_ref, var_ref))
  
  # fit measurements
  err_tp=min(dtp[dtp != 0], na.rm=TRUE)/10.
  msp=bspline::fitsmbsp(tp, mf[, -1L, drop=FALSE], n=nsp, nki=nki, monotone=mono, positive=1, lieq=lieq, control=list(monotone=TRUE, errx=err_tp, trace=0), estSD=TRUE)
  # remove metab's NA
  ina=names(which(apply(bspline::bsppar(msp)$qw, 2, function(vc) anyNA(vc))))
  if (length(ina)) {
    #browser()
    sto=stofull[-match(ina, rownames(stofull)),,drop=FALSE]
    mf=mf[,-match(ina, colnames(mf)),drop=FALSE]
    e=environment(msp)
    ibad=match(ina, colnames(e$qw))
    e$qw=e$qw[,-ibad,drop=FALSE]
    e$sdy=e$sdy[-ibad]
    e$sdqw=e$sdqw[,-ibad,drop=FALSE]
  } else {
    sto=stofull
  }
  parm=bspline::bsppar(msp)
  nmet=nrow(sto)
  nmetfull=nrow(stofull)
  nflux=ncol(sto)
  nwm=nrow(parm$qw)
  # first derivatives
  dsp=bspline::dbsp(msp)
  pard=bspline::bsppar(dsp)
  qwd=pard$qw
  nwv=nrow(pard$qw)
  # complete derivatives by 0 for non measured metabs
  qwd0=matrix(0., nrow(qwd), nmet)
  colnames(qwd0)=rownames(sto)
  qwd0[,colnames(qwd)]=qwd

  # calculates weights for fluxes
  if (!dls) {
    # build integral LS
    # qwm0: metab coeffs including 0s, it will be rhs
    qwm0=matrix(0., nrow=nwm, ncol=nmet)
    colnames(qwm0)=rownames(sto)
    qwm0[,colnames(parm$qw)]=parm$qw
    # unknowns: qwv and mstart
    mstart=setNames(double(nmet), colnames(qwm0))
    mstart[colnames(parm$qw)]=parm$qw[1,] # starting values where known
    # integration matrix (cumsum(...))
    ima=diag(nrow=nwm-1)
    ima=(row(ima)>=col(ima))+0.
    ima=rbind(0., arrApply::arrApply(ima, 2, "multv", v=diff(pard$xk, lag=nsp)/nsp))
    # f2m: ima%*%qwv%*%t(sto)+rep(1,nkmet)%o%mstart
    jqw=aperm(ima%o%t(sto), c(1L,4L,2L,3L))
    dim(jqw)=c(prod(dim(jqw)[1L:2L]), prod(dim(jqw)[3L:4L]))
    jadd=rep(1., nrow(ima))%o%diag(nrow=nmet)
    dim(jadd)=c(nrow(ima)*nmet, nmet)
    jtot=cbind(jadd, jqw)
    # inequalities: qwm >= 0; dm/dt monotone for some metabs
    if (any(monotone != 0)) {
      if (length(monotone) == 1L) {
        mono=setNames(rep(monotone, nmet), rownames(sto))
      } else {
        mono=monotone[rownames(sto)]
      }
      i=names(which(mono > 0))
      umono=aperm(diag(nrow=nwv)%o%sto[i,,drop=FALSE], c(1L,3L,2L,4L))
      dim(umono)=c(nwv*length(i), nwv*nflux)
      i=names(which(mono < 0))
      tmp=aperm(diag(nrow=nwv)%o%(-sto[i,,drop=FALSE]), c(1L,3L,2L,4L))
      dim(tmp)=c(nwv*length(i), nwv*nflux)
      umono=cbind(matrix(0., nrow=nrow(umono)+nrow(tmp), ncol=nmet), rbind(umono, tmp))
      cmono=double(nrow(umono))
    } else {
      umono=NULL
      cmono=NULL
    }
    # positivity constraint matrix u is the same as jtot
    u=rbind(jtot, umono)
    co=double(nrow(u))
    # solve ILS
    st=system.time({p=nlsic::lsi(jtot, c(qwm0), u=u, co=co)})
    # extract mst (metabolite starting values) and qwv
    icnst=seq_len(nmet)
    mst=p[icnst]
    names(mst)=rownames(sto)
    qwv=p[-icnst]
    dim(qwv)=c(nwv, nflux)
    colnames(qwv)=colnames(sto)
    # find sd
    ## cov of rhs
    covm0=slam::simple_triplet_zero_matrix(prod(dim(qwm0)))
    nqr=nrow(qwm0)
    iqr=seq_len(nqr)
    for (m in colnames(parm$qw)) {
      im=match(m, colnames(qwm0))
      i=(im-1L)*nqr+iqr
      covm0[i,i]=parm$covqw*(parm$sdy[m]**2)
    }
    ## sdv of jtot
    s=base::svd(jtot) # it is supposed to be full rank, otherwise lsi() would stop already
    jinv=tcrossprod(arrApply::arrApply(s$v, 2, "multv", v=1./s$d), s$u) # generalized inverse of jtot
    sdp=sqrt(diag(tcrossprod(slam::matprod_simple_triplet_matrix(jinv, covm0),jinv)))
    # extract sd of mst and qwv
    sdmst=sdp[icnst]
    sdfl=sdp[-icnst]
    dim(sdfl)=c(nwv, nflux)
    colnames(sdfl)=colnames(sto)
  } else {
    # differential LS
    # generalized sto inverse
    st=system.time({
    s=svd(sto)
    d=s$d
    d[d <= d[1L]*1.e-10]=0.
    d[d != 0.]=1./d[d != 0.]
    stoinv=s$v%*%(d*t(s$u))
    dimnames(stoinv)=rev(dimnames(sto))
    # solve dLS
    qwv=tcrossprod(qwd0, stoinv)
    })
    # find sd
    ## estimate cov of qwd
    dm=bspline::dmat(nrow(parm$qw), parm$xk, parm$n) # diff matrix
    sdd=sqrt(diag(tcrossprod(dm%*%parm$covqw, dm)))%o%parm$sdy
    #sdd=sdd[,!is.na(parm$sdy),drop=FALSE]
    sdd0=matrix(0., nrow=nrow(sdd), ncol=nmet)
    colnames(sdd0)=rownames(sto)
    sdd0[,colnames(sdd)]=sdd
    sit=t(stoinv)
    sdfl=matrix(0., nrow=nrow(sdd0), ncol=nflux)
    for (i in seq_len(nrow(sdd))) {
      di=sdd0[i,]
      dsit=di*sit
      for (j in seq_len(ncol(sto)))
        sdfl[i,j]=sqrt(dsit[,j]%*%dsit[,j])
    }
  }
  #print(st)
  # rate splines
  vsp=bspline::par2bsp(nsp-1L, qwv, pard$xk)
  e=environment(vsp)
  e$sdqw=sdfl # store sd
  # metab restored from S*v
  qwde=qwv%*%t(stofull) # dm/dt estimated from rates
  fsp=bspline::par2bsp(nsp-1L, qwde, pard$xk)
  const=setNames(double(nmetfull), colnames(qwde))
  if (!dls) {
    const[names(mst)]=mst
    #browser()
  } else {
    const[colnames(parm$qw)]=parm$qw[1,] # starting values where known
  }
  isp=bspline::ibsp(fsp, const=const)
  pari=bspline::bsppar(isp)
  if ((dls || length(ina) > 0) && any({imin <- apply(pari$qw, 2, min); ineg <- imin < 0.})) {
    e=environment(isp)
    #browser()
    ineg=names(ineg[ineg])
    e$qw[,ineg]=arrApply::arrApply(e$qw[,ineg,drop=FALSE], 2, "addv", v=-imin[ineg])
  }
  # residuals dm/dt-sto*f
  rsp=bspline::par2bsp(nsp-1L, qwd0-qwde[,colnames(qwd0)], pard$xk)
  # integrated residuals m-\int sto*f
  qw0=matrix(0., nwm, nmetfull)
  colnames(qw0)=colnames(pari$qw)
  qw0[,colnames(parm$qw)]=parm$qw
  risp=bspline::par2bsp(nsp, qw0-pari$qw, parm$xk)
  #browser()
  # chi2 test
  rss=colSums(as.matrix(mf[,-1L]-isp(mf$Time, colnames(mf)[-1L]))**2, na.rm=TRUE)
  chi2=rss/var_ref
  df=colSums(!is.na(mf[,-1L]))-(if(dls) nrow(qwv) else nwm)
  pval=stats::pchisq(chi2, df=df, lower=FALSE)
  chi2tab=data.frame(rss=rss, var_ref=var_ref, df=df, chi2=chi2, pval=pval)
  #browser()
  
  res=list(mf=mf, tp=tp, tpp=tpp, sto=sto, invsto=if (dls) stoinv else NULL, stofull=stofull, msp=msp, vsp=vsp, fsp=fsp, dsp=dsp, isp=isp, rsp=rsp, risp=risp, sdrate=sdfl, chi2tab=chi2tab)
  # atom balance
  if (length(atomlen)) {
    asp=bspline::par2bsp(nsp, colSums(t(parm$qw)*atomlen[colnames(parm$qw)]), parm$xk)
    pari=bspline::bsppar(isp)
    iasp=bspline::par2bsp(nsp, colSums(t(pari$qw)*atomlen[colnames(pari$qw)]), pari$xk)
    res$asp=asp
    res$iasp=iasp
    res$atomlen=atomlen
  }
  res
}
