#' square of L2 norm
#' @noRd
#' @keywords internal
norm2=function(x) sum(x*x, na.rm=TRUE)

#' matrix row-vector multiplication term-by-term.
#' @noRd
#' @keywords internal
mrowv=function(m,v) arrApply::arrApply(m, 2, "multv", v=v)

#' Generalized inverse matrix based on svd
#' @noRd
#' @keywords internal
sinv=function(a, s=base::svd(a), tol=1.e-10) {
  d=s$d
  d[d <= d[1L]*tol]=0.
  d[d != 0.]=1./d[d != 0.]
  inv=tcrossprod(mrowv(s$v, d), s$u) # generalized inverse
  dimnames(inv)=rev(dimnames(a))
  inv
}
#' Generalized inverse matrix based on svd for sd calculation
#' @noRd
#' @keywords internal
sinvsd=function(a, s=base::svd(a), tol=1.e-10) {
  d=s$d
  d[]=ifelse(d <= d[1L]*tol, 1./(d[1L]*tol), 1./d)
  inv=tcrossprod(mrowv(s$v, d), s$u) # generalized inverse
  dimnames(inv)=rev(dimnames(a))
  inv
}
#' Generalized inverse of full-rank matrix (based on QR)
#' @noRd
#' @keywords internal
qrinv=function(a, qa=base::qr(a, LAPACK=TRUE), tol=1.e-10) {
  d=abs(diag(qa$qr))
  irank=seq_len(sum(d >= d[1L]*tol))
  backsolve(qr.R(qa)[irank, irank], t(qr.Q(qa)[,irank]))[qa$pivot[irank],,drop=FALSE]
}
#' Replace abs <= tol by 0
#' @noRd
#' @keywords internal
tol0=function(x, tol=1.e-10) {x[abs(x) <= tol]=0; x}
#' Replace NA by 0
#' @noRd
#' @keywords internal
na.0=function(x) {x[is.na(x)]=0; x}
#' Drop abs <= tol
#' @noRd
#' @keywords internal
dr0=function(x, tol=1.e-10) {x[abs(x) >= tol]}
#' Iterate on name and item in named list
#' @noRd
#' @keywords internal
nm_lapply=function(X, FUN, ...) {
  nms=names(X)
  stopifnot(length(nms) == length(X))
  setNames(lapply(seq_along(X), function(i) FUN(nms[i], X[[i]], ...)), nms)
}
#' sdy for column covariance of \code{qw%*%A}
#' @noRd
#' @keywords internal
sdyA=function(sdy, A, transp=FALSE) {
  # sdyA =sqrt( diag (A' diag(sdy**2) A)) if transp is FALSE
  if (transp)
    sqrt(diag(tcrossprod(mrowv(A, sdy^2), A)))
  else
    sqrt(diag(crossprod(A, (sdy^2)*A)))
}

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
#'   matpoints(tp, res$mf[,-1], pch=".", cex=0.5)
#'   legend("topright", legend=colnames(bsppar(res$msp)$qw), lty=1:5, col=1:6, cex=0.75)
#'   # plot rates
#'   dev.new()
#'   matplot(tpp, res$vsp(tpp), type="l")
#'   ref=t(read.delim(file.path(ddir, "glyco_teusink.flux.tsv"), row.names=1, check.names=FALSE))
#'   tf=as.numeric(rownames(ref)) # reference rate time points
#'   nm_rate=colnames(bsppar(res$vsp)$qw)
#'   itf=(tf >= min(tp) & tf <= max(tp))
#'   matpoints(tf[itf], ref[itf, nm_rate, drop=FALSE], pch=".", cex=0.5)
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
#' @importFrom stats setNames na.omit
#' @importFrom graphics legend matlines matpoints polygon points lines
#' @importFrom utils packageVersion
#' 
#' @export
cli=function(args=commandArgs(trailingOnly=TRUE)) {
  Specie <- Time <- Value <- NULL # to keep 'R CMD check' calm about subset() arguments
  if (TRUE) {
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
        "Internal knot number for B-splines fitting specie and rate
        dynamics, must be less or equal to length(time_points)-2
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
    realistic fit if corresponding 'Value' is set to -1. Valid only for ILS. Cf. also options --increasing and --decreasing."
      ),
      make_option(c("--increasing"), type="character", default="",
        help="List of coma separated species that are supposed to be monotonously increasing, e.g. '--increasing=LAC,ETOH'. If a specie is present both in file MONO and in this option,
    this option takes the precedence. Valid only for ILS."
      ),
      make_option(c("--decreasing"), type="character", default="",
        help="List of coma separated species that are supposed to be monotonously decreasing, e.g. '--decreasing=GLC'. If a specie is present both in file MONO and in this option,
    this option takes the precedence. Valid only for ILS."
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
        help="use Differential Least Squares formulation (default: FALSE, i.e. Integral Least Squares (ILS) are used.)"
      ),
      make_option(c("--wsd"), action="store_true", default=FALSE,
        help="weight Least Squares by square root of covariance matrix (default: FALSE). Not compatible with under-detremined stoichiometric matrix, i.e. when there are no sufficient specie dynamics measurements."
      ),
      make_option(c("--npi"), type="integer", default=300L, help=
        "Number of plot intervals for smooth curve plotting
    [default %default]"
      ),
      make_option(c("--sf"), type="character", default="", help=
        "Names of species separated by comma ',' for which scailng factors must be estimated.
    [default %default]"
      ),
      make_option(c("--nosf"), type="character", default="", help=
        "Names of species separated by comma ',' for which scailng factors must NOT be estimated but they are estimated for all others. This option cannot be used simultaneously with '--sf'.
    [default %default]"
      ),
      make_option(c("--sderr"), type="character", default="", help=
        "Couples of species and experimental error SD separated by comma ',' e.g. 'GLC=0.01,F6P=0.02'"
      ),
      make_option(c("--fsd"), type="double", default=2., help=
        "SD factor for plotting gray band \u00b1fsd*SD around spline curves. Use '--fsd=0' to cancel these bands. [default %default]"
      ),
      make_option(c("--regular_grid"), action="store_true", default=TRUE, help=
        "use regular knot grid (default: TRUE)]"
      ),
      make_option(c("--pch"), type="character", default=".", help=
        "Plotting character for dots. [default %default]"
      )
    )
  }
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
  ikeep=rep(TRUE, nrow(mf))
  if (opt$skip > 0L)
    ikeep[seq_len(opt$skip)]=FALSE
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
  v=strsplit(opt$increasing, ",")[[1L]]
  if (any(ibad <- !(v %in% names(mono)))) {
    vclose=sapply(v[ibad], function(val) paste0(val, " (did you mean ", paste0(agrep(val, x=names(mono), value=TRUE), collapse=", "), ")"))
    stop("The following metabolites asked to be increasing, are not found in the network:\n\t", paste0(vclose, collapse="\n\t"))
  }
  mono[v]=1
  v=strsplit(opt$decreasing, ",")[[1L]]
  if (any(ibad <- !(v %in% names(mono)))) {
    vclose=sapply(v[ibad], function(val) paste0(val, " (did you mean ", paste0(agrep(val, x=names(mono), value=TRUE), collapse=", "), ")"))
    stop("The following metabolites asked to be decreasing, are not found in the network:\n\t", paste0(vclose, collapse="\n\t"))
  }
  mono[v]=-1
  # read atom lengths
  if (!is.null(opt$atom)) {
    datom=read.delim(opt$atom, comment.char="#")
    atomlen=setNames(double(nrow(sto)), rownames(sto))
    atomlen[datom[,1L]]=datom[,2L]
  } else {
    atomlen=NULL
  }
  # prepare NA list from --lna
  lna=unique(sapply(strsplit(opt$lna, ",")[[1L]], trimws))
  if (any(ibad <- !(lna %in% rownames(sto)))) {
    vclose=sapply(lna[ibad], function(val) paste0(val, " (did you mean ", paste0(agrep(val, x=rownames(sto), value=TRUE), collapse=", "), "?)"))
    stop("Following species were labeled as NA but are not present in the network:\n\t", paste0(vclose, collapse="\n\t"))
  }
  ina=colnames(mf) %in% lna
 
  # prepare nmsf vector from --sf or --nosf
  nmsf=character(0L)
  if (nchar(opt$sf)) {
    if (nchar(opt$nosf))
      stop("Both of '--sf' (", opt$sf, ") and '--nosf' (", opt$nosf, ") cannot be used simultaneously.")
    nmsf=unique(sapply(strsplit(opt$sf, ",")[[1L]], trimws))
    if (any(ibad <- !(nmsf %in% rownames(sto)))) {
      vclose=sapply(nmsf[ibad], function(val) paste0(val, " (did you mean ", paste0(agrep(val, x=rownames(sto), value=TRUE), collapse=", "), "?)"))
      stop("Following species were asked for scaling factors but are not present in the network:\n\t", paste0(vclose, collapse="\n\t"))
    }
  }
  if (nchar(opt$nosf)) {
    nmsf=unique(sapply(strsplit(opt$nosf, ",")[[1L]], trimws))
    if (any(ibad <- !(nmsf %in% rownames(sto)))) {
      vclose=sapply(nmsf[ibad], function(val) paste0(val, " (did you mean ", paste0(agrep(val, x=rownames(sto), value=TRUE), collapse=", "), "?)"))
      stop("Following species were asked for NOT be used with scaling factors but are not present in the network:\n\t", paste0(vclose, collapse="\n\t"))
    }
    nmsf=setdiff(rownames(sto), nmsf)
  }
  nmsf=setdiff(nmsf, lna)
  if (all(rownames(sto) %in% c(nmsf, lna)))
    stop("Scaling factors are requested for all available species which is meaningless. At least one specie must be excluded from --sf list")
  #print(c("opt=", opt))
  
  # prepare sderr
  sderr=double(0L)
  if (nchar(opt$sderr)) {
    sdtmp=sapply(strsplit(opt$sderr, ",")[[1L]], trimws)
    sdtmp=lapply(strsplit(sdtmp, "="), trimws)
    if (any(lengths(sdtmp) != 2L))
      stop("Badly formatted --sdref '", opt$sderr, "'. Expecting couples 'SPECIE1=VALUE1,...,SPECIEN=VALUEN'")
    sderr=setNames(sapply(sdtmp, `[[`, 2L), sapply(sdtmp, `[[`, 1L))
    suppressWarnings(storage.mode(sderr) <- "double")
    if (anyNA(sderr))
      stop("A value could not be converted to real number in --sdref '", opt$sderr, "'. Expecting couples 'SPECIE1=VALUE1,...,SPECIEN=VALUEN'")
    if (any(ibad <- !(names(sderr) %in% rownames(sto)))) {
      vclose=sapply(names(sderr)[ibad], function(val) paste0(val, " (did you mean ", paste0(agrep(val, x=rownames(sto), value=TRUE), collapse=", "), "?)"))
      stop("Following species were assigned error SD in '", opt$sderr, "' but are not present in the network:\n\t", paste0(vclose, collapse="\n\t"))
    }
    if (any(names(sderr) %in% lna))
      sderr=sderr[!names(sderr) %in% lna]
  }
  
  # main call
  res=fdyn(mf[ikeep,!ina], sto, nsp=opt$norder, nki=opt$knot, lieq=lieq[setdiff(names(lieq), lna)], monotone=mono, dls=opt$dls, atomlen=atomlen, npi=opt$npi, wsd=opt$wsd, nmsf=nmsf, sderr=sderr, regular_grid=opt$regular_grid)
  res$mffull=mf
  res$opt=opt
  res$date=date()
  #res
  # write result files (rd is a temporary dir for results)
  # at the end we'll move all files into a zip archive in the working dir
  #rd=tempfile(pattern="fdyn")
  if (is.null(opt$out)) {
    rd=tools::file_path_sans_ext(opt$meas)
  } else {
    rd=opt$out
  }
  if (nchar(rd)) { # write results
    if (!dir.exists(rd))
      dir.create(rd)
    tp=res$tp
    tpp=res$tpp
    ratp=range(tp, na.rm=TRUE)
    tppr=c(tpp, rev(tpp)) # for SD band plotting
    bcol=do.call(rgb, as.list(c(col2rgb(1)/255, 0.3))) # band color
    # define plot function
    plotsp=function(fname, sp, main, ylab, data=NULL, sto=NULL, pch=opt$pch) {
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
        matpoints(tp, data, pch=pch, cex=0.5)
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
          points(tp, d, pch=pch, cex=0.5)
        lines(tpp, mc[,m], lwd=1.5)
        if (!is.null(p$sdqw))
          polygon(tppr, c(msdp[,m], rev(msdm[,m])), border=NA, col=bcol)
        if (!is.null(sto)) {
          matlines(tpp, fl, lty=c(2:5,1), col=c(2:6,1), lwd=1.5)
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
      datom=colSums(t(res$mf[,-1L])*atomlen[colnames(res$mf)[-1L]], na.rm=TRUE)
      ac=res$asp(tpp)
      iac=res$iasp(tpp)
      plot(1, xlim=ratp, main="Atom balance evolution", ylim=range(ac, iac, datom, na.rm=TRUE),
        xlab="Time", ylab="Total atom number", t="n", lwd=1.5)
      matlines(tpp, cbind(ac, iac), lwd=1.5)
      points(tp, datom, pch=opt$pch, cex=0.5)
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
#browser()
    grDevices::cairo_pdf(file.path(rd, "imet%03d.pdf"))
    inames=colnames(bspline::bsppar(res$isp)$qw)
    tmp=matrix(NA_real_, nrow(mf[ikeep,,drop=FALSE]), length(inames)) # inject here mf values
    colnames(tmp)=inames
    cnm=intersect(inames, colnames(res$mf)[-1L]) # NA columns in mf are absent in mc
    tmp[,cnm]=as.matrix(res$mf[,cnm])
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
    # stats
    write.table(res$chi2tab, file=file.path(rd, "stats.tsv"), sep="\t", quote=FALSE, col.names=NA)
    # Readme.md
#browser()
    if (length(nmsf)) {
      cat(file=file.path(rd, "sf.tsv"), sep="",
        "Specie\tValue\n")
      write.table(as.matrix(res$sf), file=file.path(rd, "sf.tsv"), append=TRUE, sep="\t", quote=FALSE, col.names=FALSE)
    } else {
      unlink(file.path(rd, "sf.tsv"))
    }
    
    cat(file=file.path(rd, "Readme.md"), sep="\n",
      "# Retrieving reaction rate dynamics from specie kinetics (dynafluxr results)", "",
      paste0("This is the result files produced by dynafluxr R package (v", utils::packageVersion("dynafluxr"), ") on ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z (%Z).")), "",
      "The command to reproduce these results is:", "",
      paste0("`Rscript --vanilla -e 'dynafluxr::cli()' ", paste0(shQuote(args), collapse=" "), "`"),
      "", "## File contents", "",
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
      " - `stats.tsv`: table with chi2 tests per compound;",
      paste0(" - `env.RData`: stored R list `res` such as returned by `dynafluxr::fdyn()`. It can be read in R session with `e=new.env(); load('", file.path(rd, "env.RData"), "', envir=e)` and then used to retrieve e.g. integrated compounds as `icmpnd=e$res$isp(e$res$tpp)`;"),
      " - `Readme.md`: this file;",
      if (length(nmsf)) " - `sf.tsv`: estimated scaling factors;" else ""
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
#' @param npi Integer scalar, indicates a number of plot intervals to produce smooth plots.
#'   (default 300)
#' @param wsd Logical scalar, if TRUE, indicates that differential least squares
#'   should be resolved with residuals weighted by a factor of covariance matrix.
#'   (default FALSE, i.e. no weighting is used)
#' @param nmsf Character vector, list of species for which scaling factor maust be estimated for --dls.
#' @param tol Double scalar, tolerance for detecting singular matrices and solving linear systems
#' @param regular_grid Logical scalar, use regular knot grid (default: TRUE)
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
#'   \item{stofull:}{ stoichiometric matrix before a possible NA elimination}
#'   \item{stoinv:}{ pseudo-inverse of sto}
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
#'   \item{sf:}{ named scale factor vector}
#'   \item{internal_knot_ref:}{ number of internal knots used for estimation of var_ref}
#' }
#' @importFrom bspline fitsmbsp dbsp bsppar par2bsp ibsp
#' @importFrom slam simple_triplet_zero_matrix matprod_simple_triplet_matrix
#' @importFrom arrApply arrApply
#' @importFrom stats var pchisq
#' @importFrom nlsic lsi
#' @export
fdyn=function(mf, stofull, nsp=4L, nki=5L, lieq=NULL, monotone=0, dls=FALSE,
    atomlen=NULL, npi=300L, wsd=FALSE, nmsf=character(0L), sderr=NULL, tol=1.e-10, regular_grid=TRUE) {
  tp=mf$Time
  dtp=diff(tp)
  np=length(tp)
  stopifnot(nki <= np-2L)
  tpp=seq(tp[1L], tp[np], length.out=npi+1L)
  # eliminate all NA columns from mf
  icna=apply(mf, 2L, function(v) all(is.na(v)))
  mf=mf[,!icna, drop=FALSE]
  # prepare mono_mf
  mono_mf=if (length(monotone) > 1L) monotone[colnames(mf)[-1L]] else monotone
#browser()
  # estimate var_ref with biggest d2 in variance (for chi2 test)
  control=list(errx=min(dtp[dtp != 0], na.rm=TRUE)/100., trace=0)
  nki_test=seq(max(0L, nki-5L), min(nki+5L, np-3L))
  var_test=sapply(nki_test, function(k) {
    if (!regular_grid) {
      xki=bspline::iknots(seq(tp[1L], tp[np], length.out=np),
        mf, n=nsp, nki=k, lenfit=0L, smooth=TRUE)
      s=bspline::fitsmbsp(tp, mf[, -1L, drop=FALSE], n=nsp, xki=xki, monotone=mono_mf, positive=1, control=control, lieq=lieq, estSD=FALSE, regular_grid=regular_grid)
    } else {
      s=bspline::fitsmbsp(tp, mf[, -1L, drop=FALSE], n=nsp, nki=k, monotone=mono_mf, positive=1, control=control, lieq=lieq, estSD=FALSE, regular_grid=regular_grid)
    }
    if (anyNA(bspline::bsppar(s)$qw))
      stop("NA appeared in preliminary spline fits")
    base::colMeans((s(tp)-mf[, -1L, drop=FALSE])**2L, na.rm=TRUE)
  })
  # detect maximal curvature -> ku=2*d2/(1+d1^2)^1.5
  ly=log(base::colSums(var_test, na.rm=TRUE))
  ku=2.*base::diff(ly, difference=2L)/(1.+(base::diff(ly, difference=1L, lag=2L)*0.5)^2L)^1.5
  i_ref=which.max(ku)+1L
  var_ref=var_test[,i_ref]
  var_ref[names(sderr)]=sderr**2L
#browser()

  # fit measurements
  #xki=bspline::iknots(tp, mf[, -1L, drop=FALSE], n=nsp, nki=nki, lenfit=11L)
  #xki=seq(tp[1L], tp[np], length.out=nki+2L)[c(-1L, -(nki+2L))]
  if (!regular_grid) {
    xki=bspline::iknots(seq(tp[1L], tp[np], length.out=np),
        mf, n=nsp, nki=nki, lenfit=0L, smooth=TRUE)
    msp=bspline::fitsmbsp(tp, mf[, -1L, drop=FALSE], n=nsp, xki=xki, monotone=mono_mf, positive=1, control=control, lieq=lieq, estSD=TRUE, regular_grid=regular_grid)
  } else {
    msp=bspline::fitsmbsp(tp, mf[, -1L, drop=FALSE], n=nsp, nki=nki, monotone=mono_mf, positive=1, control=control, lieq=lieq, estSD=TRUE, regular_grid=regular_grid)
  }
#browser()
  # error on metab's NA
  ina=names(which(apply(bspline::bsppar(msp)$qw, 2, function(vc) anyNA(vc))))
  if (length(ina)) {
    #browser()
    stop("Following species could not be fitted: '", paste0(ina, collapse="', '"), "'")
#    ima=match(ina, rownames(stofull))
#    ibad=which(is.na(ima))
#    if (length(ibad))
#      stop("Following names passed to --lna are not recognized: '", paste0(ina[ibad], collapse="', '"), "'")
#    sto=stofull[-ima,,drop=FALSE]
#    mf=mf[,-match(ina, colnames(mf)),drop=FALSE]
#    e=environment(msp)
#    ibad=match(ina, colnames(e$qw))
#    e$qw=e$qw[,-ibad,drop=FALSE]
#    e$sdy=e$sdy[-ibad]
#    e$sdqw=e$sdqw[,-ibad,drop=FALSE]
  }
  # keep only valid metab in sto
  sto=stofull[intersect(rownames(stofull), colnames(mf)[-1L]),,drop=FALSE]
  nrate=ncol(sto)
  nmet=nrow(sto)
  nmetfull=nrow(stofull)
  #if (nmet < nrate)
  #  stop("Number of species in stoichiometric matrix (", nrow(sto), ") is lower than reaction number (", ncol(sto), ").")
  # check rank and make stoinv
  s=svd(sto)
  d=s$d
  srank=sum(d > d[1L]*tol)
#browser()
  if (srank < nrate) {
    # minimal norm solution
    qs=qr(sto, LAPACK=TRUE)
    warning("Stoichiometric matrix rank (", srank, ") is lower than reaction number (", nrate, ").\n",
      "The solution of minimal norm will be provided.\n",
      "Undefined rates could be, for example:\n\t", paste0(colnames(sto)[tail(qs$pivot, -srank)], collapse="\n\t"))
  }
  stoinv=sinv(sto, s, tol=tol)
  stoinvsd=sinvsd(sto, s, tol=tol)
  parm=bspline::bsppar(msp)
  nwm=nrow(parm$qw)
  xk=parm$xk
  # first derivatives
  dsp=bspline::dbsp(msp)
  pard=bspline::bsppar(dsp)
  qwd=pard$qw
  nwv=nrow(pard$qw)
  # complete derivatives by 0 for non measured metabs
  qwd0=matrix(0., nrow(qwd), nmet)
  colnames(qwd0)=rownames(sto)
  qwd0[,colnames(qwd)]=qwd
  
  nb_sf=length(nmsf)
  sf=setNames(rep(1., nb_sf), nmsf)
  if (nb_sf > 0L && any(ibad <- !(nmsf %in% rownames(sto))))
    stop("following species were asked with scaling factors but are not measured:\n\t", paste0(nmsf[ibad], collapse="\n\t"))
  i_sf=match(nmsf, colnames(qwd0))
  # inequalities: qwm >= 0; dm/dt monotone for some metabs
  if (any(monotone != 0)) {
    if (length(monotone) == 1L) {
      mono=setNames(rep(monotone, nmet), rownames(stofull))
    } else {
      mono=monotone[rownames(stofull)]
    }
    i=names(which(mono > 0))
    umono=aperm(diag(nrow=nwv)%o%stofull[i,,drop=FALSE], c(1L,3L,2L,4L))
    dim(umono)=c(nwv*length(i), nwv*nrate)
    i=names(which(mono < 0))
    tmp=aperm(diag(nrow=nwv)%o%(-stofull[i,,drop=FALSE]), c(1L,3L,2L,4L))
    dim(tmp)=c(nwv*length(i), nwv*nrate)
    umono=cbind(matrix(0., nrow=nrow(umono)+nrow(tmp), ncol=nmet), rbind(umono, tmp))
    cmono=double(nrow(umono))
  } else {
    umono=NULL
    cmono=NULL
  }
  
  x=NULL
  # calculates B-splines coefs for reaction rates -> qwv
  if (!dls) {
    # build Integral LS
    # qwm0: metab coeffs including 0s, it will be rhs
    qwm0=matrix(0., nrow=nwm, ncol=nmet)
    colnames(qwm0)=rownames(sto)
    qwm0[,colnames(parm$qw)]=parm$qw
    qwm0sf=qwm0
    # unknowns: qwv and mstart
    # starting values where known
    mstart=setNames(na.0(parm$qw[1L,][colnames(qwm0)]), colnames(qwm0))
    icnst=seq_len(nmet)
    # integration matrix (cumsum(...))
    ima=diag(nrow=nwm-1)
    ima=(row(ima)>=col(ima))+0.
    ima=rbind(0., mrowv(ima, diff(pard$xk, lag=nsp)/nsp))
    # differentiate matrix
    dma=bspline::dmat(f=msp)
    qwv=tcrossprod(qwd0, stoinv)
    rownames(qwv)=paste0("k", seq_len(nrow(qwd0)))
    # rate2metab: ima%*%qwv%*%t(sto)+rep(1,nkmet)%o%mstart ~> qwm0
    jqw=aperm(ima%o%t(sto), c(1L,4L,2L,3L))
    dim(jqw)=c(prod(dim(jqw)[1L:2L]), prod(dim(jqw)[3L:4L]))
    jadd=rep(1., nrow(ima))%o%diag(nrow=nmet)
    dim(jadd)=c(nrow(ima)*nmet, nmet)
    jtot=cbind(jadd, jqw)
    resj=qwm0
    # positivity constraint matrix u is the same as jtot
    u=rbind(jtot, umono)
    co=double(nrow(u))
#browser()
    # prepare chol factor of cov matrix (t(R)%*%R) for weighting
    if (wsd) {
      if (srank < nrate)
        stop("Option '--wsd' is not compatible with under-determined stoichiometric matrix")
      rcovqw0=chol(parm$covqw)
      rcovqw0=t(backsolve(rcovqw0, diag(ncol(rcovqw0))))
      isdy0=structure(rep(1./mean(parm$sdy), nmet), names=rownames(sto))
      isdy0[names(parm$sdy)]=1./parm$sdy;
      # weight jtot and qwm0 by rcov
      for (im in seq_len(nmet)) {
        ir=seq(nwm)+nwm*(im-1L)
        jtot[ir,]=isdy0[im]*(rcovqw0%*%jtot[ir,])
        resj[,im]=isdy0[im]*(rcovqw0%*%resj[,im])
      }
    }
    if (nb_sf) {
      # prepend scaling factor to unknown vector
      jsf=array(0., dim=c(dim(qwm0), nb_sf))
      for (i in seq_len(nb_sf))
        jsf[,i_sf[i],i]=-resj[,nmsf[i]]
      jsf=matrix(jsf, prod(dim(resj)), nb_sf)
      jtot0=jtot # to remove
      jtot=cbind(jsf, jtot)
      resj[,nmsf]=0.
      u0=u
      u=cbind(matrix(0., nrow(u), nb_sf), u)
    }
    # solve ILS
    colnames(jtot)=c(paste0("sf.", nmsf, recycle0=TRUE),
      paste0("cnst.", rownames(sto)),
      t(outer(colnames(qwv), seq_len(nwv), paste, sep=".")))
#browser()
    # build equality constraints if asked: matrix -> mate, rhs vector -> ce
    mate=NULL
    ce=NULL
    if (length(lieq)) {
      # mate: neq x (nmet + (nrow*ncol)_qwv)
      # we need a path from qwm0 to qwv which is the same as in jtot
      # e%*%qwm0[,met]=ce => e%*%jtot[imet,]%*%(const,qwv)=ce =>
      # mete = e%*%jtot[imet,]
      # for a given equality: $e,$ce are a matrix and rhs vector for qwm0[,met]
      nqr=nrow(qwm0)
      iqr=seq_len(nqr)
      liece=nm_lapply(lieq, function(m, meq) {
          if (length(meq) == 0L) {
              NULL
          } else {
              im=match(m, colnames(qwm0))
              i=(im-1L)*nqr+iqr
              list(e=bspline::bsc(meq[,1L], n=nsp, xk)%*%jtot[i,,drop=FALSE], ce=meq[,2L])
          }
      })
      mate=Reduce(rbind, lapply(liece, `[[`, "e"), NULL)
      ce=Reduce(c, lapply(liece, `[[`, "ce"), NULL)
      st=system.time({p=nlsic::lsie_ln(jtot, c(resj), u=u, co=co, e=mate, ce=ce, rcond=1./tol)})
    } else {
      st=system.time({p=nlsic::lsi_ln(jtot, c(resj), u=u, co=co, rcond=1./tol)})
    }
    if (anyNA(p))
      stop("Error in least squares with constraints:\n", attr(p, "mes"))
    if (nb_sf) {
      sf[]=p[seq_len(nb_sf)]
      p=p[-seq_len(nb_sf)]
    }
    # extract mst (metabolite starting values) and qwv
    mst=p[icnst]
    names(mst)=rownames(sto)
    qwv[]=p[-icnst]
    dim(qwv)=c(nwv, nrate)
    colnames(qwv)=colnames(sto)
    # find sd
    ## cov of rhs
    covm0=slam::simple_triplet_zero_matrix(prod(dim(qwm0)))
    if (wsd) {
      i=seq_len(ncol(covm0))
      covm0[cbind(i,i)]=1.
    } else {
      nqr=nrow(qwm0)
      iqr=seq_len(nqr)
      for (m in colnames(parm$qw)) {
        im=match(m, colnames(qwm0))
        i=(im-1L)*nqr+iqr
        covm0[i,i]=parm$covqw*(parm$sdy[m]**2)
      }
    }
    jinv=sinvsd(jtot) # generalized inverse of jtot for sd
#browser()
    covp=tcrossprod(slam::matprod_simple_triplet_matrix(jinv, covm0),jinv)
    sdp=sqrt(diag(covp))
    if (nb_sf) {
      i_sf=seq_len(nb_sf)
      sdp=sdp[-i_sf]
      covp=covp[-i_sf, -i_sf, drop=FALSE]
    }
    # extract sd of mst and qwv
    sdmst=sdp[icnst]
    sdconst=sdmst
    sdrate=sdp[-icnst]
    covqwv=covp[-icnst,-icnst]
    # covqwd is mean of nrate submatrices
    iwv=seq_len(nwv)
    lcov=lapply(seq(nrate)-1L, function(i) covqwv[i*nwv+iwv,i*nwv+iwv])
    covqwd=Reduce("+", lcov)/nwv
    sdyd=sqrt(sapply(lcov, function(mat) mean(mat/covqwd)))
    dim(sdrate)=c(nwv, nrate)
    colnames(sdrate)=colnames(sto)
  } else {
    # differential LS
    # generalized sto inverse
    # prepare chol factor of cov matrix (t(R)%*%R) for weighting
    if (!is.null(umono))
      stop("Monotonicity requirement is not compatible with --dls option")
    ## estimate cov of qwd
    dm=bspline::dmat(nrow(parm$qw), parm$xk, parm$n) # diff matrix
    covqwd=tcrossprod(dm%*%parm$covqw, dm)
    st=system.time({
    x0=c(rep(1., nb_sf), tcrossprod(qwd0, stoinv))
    f_rbax=function(x, ...) {
      mc=match.call()
      with(list(...), {
        nrc=lengths(dimnm)
        qwv=structure(if (nb_sf) x[-isf] else x, dim=nrc, dimnames=dimnm)
        res=if (!is.null(rit)) rit%*%tcrossprod(qwv, sto) else tcrossprod(qwv, sto)
        qwd0sf[,nmsf]=mrowv(qwd0[,nmsf, drop=FALSE], x[isf]) # qwd0sd
        if (mc[[1L]] == "f_resid") {
          res=qwd0sf-res
        } else {
          res[,nmsf]=res[,nmsf]-qwd0sf[,nmsf]
        }
        # preconditioning
        c(if (nb_sf) diag(pinv%*%res[,nmsf,drop=FALSE]) else NULL, tcrossprod(res, stoinv)) #stosdinv))
      })
    }
    if (wsd) {
      if (srank < nrate)
        stop("Option '--wsd' is not compatible with under-determined stoichiometric matrix")
      rcovqwd0=chol(covqwd)
      rcovqwd0i=backsolve(rcovqwd0, diag(ncol(rcovqwd0)), transpose=TRUE)
      isdy0=setNames(rep(1./mean(parm$sdy), nmet), rownames(sto))
      isdy0[names(parm$sdy)]=1./parm$sdy;
      qwd0sd=structure(mrowv(rcovqwd0i%*%qwd0, isdy0), dimnames=list(paste0("k", seq_len(nrow(qwd0))), names(isdy0)))
      qwv=matrix(0., nwv, ncol(sto))
      dimnames(qwv)=list(paste0("k", seq_len(nrow(qwv))), colnames(sto))
      # gmresls
      pinv=if (nb_sf) t(qr.Q(qr((qwd0sd[,nmsf,drop=FALSE])))) else matrix(0, 0L, 0L)
      qwd0sf=qwd0sd
      # solve for x=(sf, qwv)
      stosd=sto*isdy0
      qrs=qr(stosd)
      qt=t(qr.Q(qrs))
      x=gmresls::gmresls(f_rbax, f_rbax, x0=x0, nmsf=nmsf, nb_sf=nb_sf, isf=seq_len(nb_sf), dimnm=dimnames(qwv), pinv=pinv, qwd0=qwd0sd, qwd0sf=qwd0sf, sto=stosd, stoinv=qt, rit=rcovqwd0i, r=rcovqwd0)
      sf[]=x[seq_len(nb_sf)]
      qwv[]=if (nb_sf) x[-seq_len(nb_sf)] else x

      # find sd
      stosdinv=backsolve(qr.R(qrs), qt)
      sdyd=sdyA(rep(1., nrow(sto)), stosdinv, transp=TRUE)
      sdrate=structure(sqrt(diag(covqwd))%o%sdyd, dimnames=list(NULL, colnames(sto)))
    } else {
      # solve plain dLS
      #qwv=tcrossprod(qwd0, stoinv)
      qwv=structure(double(nrow(qwd0)*ncol(sto)), dim=c(nrow(qwd0), ncol(sto)), dimnames=list(paste0("k", seq_len(nrow(qwd0))), colnames(sto)))
      # estimate scaling factors (sf) if asked
      pinv=if (nb_sf) sinv(qwd0[,nmsf,drop=FALSE]) else matrix(0., 0L, 0L)
      i_sf=match(nmsf, colnames(qwd0))
      qwd0sf=qwd0
      qwd0sf[,i_sf]=mrowv(qwd0[,i_sf,drop=FALSE], sf)
      # solve for x=(sf, qwv)
      x=gmresls::gmresls(f_rbax, f_rbax, x0=x0, nmsf=nmsf, nb_sf=nb_sf, isf=seq_len(nb_sf), dimnm=dimnames(qwv), pinv=pinv, qwd0=qwd0, qwd0sf=qwd0sf, sto=sto, stoinv=stoinv, rit=NULL)
      sf[]=x[seq_len(nb_sf)]
      qwv[]=if (nb_sf) x[-seq_len(nb_sf)] else x

      #print(c("sf=", sf, "; rss=", norm2(t(qwd0)-sto%*%t(qwv))))
      # find sd
      sdyd=sdyA(parm$sdy, stoinvsd, transp=TRUE)
      sdrate=structure(sqrt(diag(covqwd))%o%sdyd, dimnames=list(NULL, colnames(sto)))
    }
    sdconst=double(nrow(sto))
    }) # system.time()
  }
  names(sdconst)=rownames(sto)
#browser()
  if (nb_sf) {
    qwd0[,i_sf]=mrowv(qwd0[,i_sf,drop=FALSE], sf)
    i=nmsf %in% colnames(parm$qw)
    environment(msp)$qw[,nmsf[i]]=mrowv(parm$qw[,nmsf[i],drop=FALSE], sf[i])
    parm=bsppar(msp)
    mf[,nmsf[i]]=mrowv(as.matrix(mf[,nmsf[i],drop=FALSE]), sf[i])
    var_ref[nmsf]=var_ref[nmsf]*(sf^2L)
  }
  #print(st)
  # rate splines
  vsp=bspline::par2bsp(nsp-1L, qwv, pard$xk)
  e=environment(vsp)
  e$sdqw=sdrate # store sd
  # metab restored from S*v
  qwde=qwv%*%t(stofull) # dm/dt estimated from rates
  fsp=bspline::par2bsp(nsp-1L, qwde, pard$xk)
  const=setNames(double(nmetfull), colnames(qwde))
  sdconstfull=const
  sdconstfull[names(sdconst)]=sdconst
#browser()
  isp=bspline::ibsp(fsp) # add const later
  if (!dls) {
    const[names(mst)]=mst
    #browser()
  } else {
    #const[colnames(parm$qw)]=parm$qw[1,] # starting values where known
    #least squares
    est=colMeans(isp(tp, colnames(mf)[-1L]))
    meas=colMeans(mf[,-1L], na.rm=TRUE)
    const[colnames(mf)[-1L]]=meas-est
  }
  ei=environment(isp)
  ei$qw=arrApply::arrApply(ei$qw, 2L, "addv", v=const)
#browser()
  if (any({imin <- apply(ei$qw, 2L, min); ineg <- imin < 0.})) {
    ineg=names(ineg[ineg])
    ei$qw[,ineg]=arrApply::arrApply(ei$qw[,ineg,drop=FALSE], 2L, "addv", v=-imin[ineg])
  }
  # sdi: SD of integrated metabolites
  # qwi = Y%*%qwde
  # SD of additive constants are ignored
  # cov(qwi) =  <Y qwde, qwde' Y'> = Y covde Y'
  Y=bspline::imat(fsp)
  covqwi=tcrossprod(Y %*% covqwd, Y)
  sdyi=sdyA(sdyd, stofull, transp=TRUE)
  sdi=sqrt(diag(covqwi)) %o% sdyi
  sdi[1L,]=sdconstfull
  ei$covqw=covqwi
  ei$sdy=sdyd
  ei$sdqw=sdi
  pari=bspline::bsppar(isp)
  
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
  rss=c(rss, Total=sum(rss))
  var_ref=c(var_ref, Total=sum(var_ref))
  chi2=rss/var_ref
  chi2[length(chi2)]=sum(head(chi2, -1L))
  df=colSums(!is.na(mf[,-1L]))-(if(dls) nrow(qwv) else nwm)
  df=c(df, Total=sum(df))
  pval=stats::pchisq(chi2, df=df, lower=FALSE)
  chi2tab=data.frame(rss=rss, var_ref=var_ref, df=df, chi2=chi2, pval=pval)
#browser()

  res=list(mf=mf, tp=tp, tpp=tpp, sto=sto, stofull=stofull, stoinv=stoinv, msp=msp, vsp=vsp,
    fsp=fsp, dsp=dsp, isp=isp, rsp=rsp, risp=risp, sdrate=sdrate, chi2tab=chi2tab,
    sf=sf, internal_knot_ref=nki_test[i_ref])
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
gui=function() {
  pkgs=c("shiny", "shinyFiles", "shinyjs")
  available=sapply(pkgs, requireNamespace, quietly=TRUE)
  if (any(!available))
    stop("Following R packages are necessary for GUI running. Install them first:\n\t", paste0(pkgs[!available], collapse="\n\t"))
  library(shiny)
  library(shinyFiles)
  library(shinyjs)
  inp2val=function(input)
    sapply(nm_par, function(nm) as.character(input[[nm]]))

  nm_par=c("knot", "norder", "lna", "increasing", "decreasing", "skip", "zip", "dls", "wsd", "npi", "sf", "nosf", "sderr", "fsd", "regular_grid", "pch")
  par_def=list(knot=5L, norder=4L, skip=0L, npi=300L, fsd=2., pch=".")
  temp_dir=tempdir()
  addResourcePath("pdfs", temp_dir)  # Make temp_dir accessible via /pdfs/
  ui <- fluidPage(
    useShinyjs(),
    titlePanel("DynaFluxR"),
    p("Choose parameters (measurement and network files are mandatory) and click on 'Run'."),
    p("Parameters can be saved in a .param file and loaded later"),
    sidebarLayout(
      sidebarPanel(
        fileInput("meas", "Measurement File (-m):"),
        fileInput("sto", "Stoichiometric Model (-s):"),
        actionButton(inputId = "btn_more", label = "more / less parameters"),
        conditionalPanel("input.btn_more % 2 == 1",
          numericInput("knot", "Knot Number (-k):", par_def$knot, min = 0),
          numericInput("norder", "Polynomial Order (-n):", par_def$norder, min = 1),
          textInput("lna", "List of NA Species (--lna):"),
          fileInput("constr", "Constraint File (-c):"),
          fileInput("atom", "Atom Length File (-a):"),
          fileInput("mono", "Monotonicity File (--mono):"),
          textInput("increasing", "Increasing Species (--increasing):"),
          textInput("decreasing", "Decreasing Species (--decreasing):"),
          numericInput("skip", "Skip First Time Points (--skip):", par_def$skip, min = 0),
          checkboxInput("zip", "Create ZIP Archive (--zip)", FALSE),
          checkboxInput("dls", "Use Differential Least Squares (--dls)", FALSE),
          checkboxInput("wsd", "Weight Least Squares (--wsd)", FALSE),
          numericInput("npi", "Plot Intervals (--npi):", par_def$npi, min = 100),
          textInput("sf", "Scaling Factors (--sf):"),
          textInput("nosf", "No Scaling Factors (--nosf):"),
          textInput("sderr", "Species Error SD (--sderr):"),
          numericInput("fsd", "SD Factor for Plot (--fsd):", par_def$fsd, min = 0),
          checkboxInput("regular_grid", "Use Regular Grid (--regular_grid)", TRUE),
          textInput("pch", "Plot Character (--pch):", par_def$pch)
        ),
        hr(),
        actionButton("run", "Run", icon=icon("bolt")),
        actionButton("stop", "Stop App"),
        hr(),
        #actionButton("save_params", "Save Parameters"),
        shinySaveButton("save_params", "Save parameters", "Choose a file to store parameters (end with .param)", multiple = FALSE),
        #fileInput("load_params", "Load Parameters (.param)"),
        shinyFilesButton("load_params", "Load Parameters (.param)", "Please select a file with parameters to load (.param)", multiple = FALSE, viewtype = "detail")
      ),
      mainPanel(
        h2("Results"),
        selectInput("pdf_select", "Select PDF:", choices = NULL),
        uiOutput("pdf_view"),
        selectInput("tsv_select", "Select TSV:", choices = NULL),
        tableOutput("tsv_table"),
        downloadButton("downloadZip", label = "Download Zip archive"),
        p(),
        p("Date: ", textOutput("date", inline=TRUE), style="font-size: x-small"),
        p("Call in Unix shell: ", style="font-size: x-small",
        span(textOutput("call", inline=TRUE), style="font-family: monospace")),
        p("R call: ", style="font-size: x-small",
          span(textOutput("callr", inline=TRUE), style="font-family: monospace")
        ),
        verbatimTextOutput("errmes"),
        div(id="legal", style="font-size: 10px;",
          tags$base(target="_blank"),
          p(
            a(id="alegal", href="javascript: void(0);", onClick="toggle('plegal')", style="text-decoration: none", target="_self", "▶"),
            strong(" Legal information about DynaFluxR")),
          p(id="plegal", style="display:none;",
            "Author: Serguei Sokol, ",
            a(href="https://www.inrae.fr/", "INRAE"),
            " / ",
            a(href="https://www.toulouse-biotechnology-institute.fr/", "TBI"),
            " / ",
            a(href="https://www.toulouse-biotechnology-institute.fr/plateformes-plateaux/cellule-mathematiques/", "Mathematics Cell"),
            br(),
            "Copyright: 2025 ",
            a(href="https://www.inrae.fr/", "INRAE"),
            "/",
            a(href="https://www.insa-toulouse.fr/", "INSA"),
            "/",
            a(href="https://www.cnrs.fr/", "CNRS"),
            br(),
            "License: ",
            a(href="https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html", "GPL-2"),
            br(),
            "Version: ", format(packageVersion("dynafluxr")),
            br(),
            "Availability: ",
            a(href="https://cran.r-project.org/package=dynafluxr", "CRAN package"),
            br(),
            "Cite paper: ", a(href="", "Later to come"),
            br(),
            "Project home: ",
            a(href="https://github.com/MetaSys-LISBP/dynafluxr", "DynaFluxR"),
            br(),
            "Issue report: ",
            a(href="https://github.com/MetaSys-LISBP/dynafluxr/issues", "github page")
          )
        ),
        tags$script('function toggle(nm_obj) {
            var obj=document.getElementById(nm_obj);
            if (obj.style.display == "block") obj.style.display = "none";
            else obj.style.display = "block";
            var alink=document.getElementById(nm_obj.replace("p", "a"));
            if (alink != undefined) {
                if (alink.innerHTML == "▼")
                    alink.innerHTML = "▶";
                else
                    alink.innerHTML = "▼";
            }
        }')
      )
    ),
    tags$script(
      "Shiny.addCustomMessageHandler('gotop', function(mes) {
           window.scrollTo(0,0);
           //console.log(mes)
        });"
    ),
    NULL
  )

  server <- function(input, output, session) {
    hd=path.expand("~")
    wd=getwd()
    volumes <- c(Home = hd, getVolumes()())
    wd=if (startsWith(wd, hd)) substring(wd, nchar(hd)+2L) else ""
    shinyFileChoose(input, "load_params", roots = volumes, session = session, defaultPath = wd)
    shinyFileSave(input, "save_params", roots = volumes, session=session, defaultPath = wd)
    observeEvent(input$run, {
      req(input$meas)
      req(input$sto)
      args <- c(
        if (!is.null(input$meas)) c("-m", input$meas$datapath),
        if (!is.null(input$sto)) c("-s", input$sto$datapath),
        if (input$knot != par_def$knot) c("-k", input$knot),
        if (input$norder != par_def$norder) c("-n", input$norder),
        if (nzchar(input$lna)) c("--lna", input$lna),
        if (!is.null(input$constr)) c("-c", input$constr$datapath),
        if (!is.null(input$atom)) c("-a", input$atom$datapath),
        if (!is.null(input$mono)) c("--mono", input$mono$datapath),
        if (nzchar(input$increasing)) c("--increasing", input$increasing),
        if (nzchar(input$decreasing)) c("--decreasing", input$decreasing),
        if (input$skip != par_def$skip) c("--skip", input$skip),
        if (input$zip == "TRUE") "--zip",
        if (input$dls == "TRUE") "--dls",
        if (input$wsd == "TRUE") "--wsd",
        if (input$npi != par_def$npi) c("--npi", input$npi),
        if (nzchar(input$sf)) c("--sf", input$sf),
        if (nzchar(input$nosf)) c("--nosf", input$nosf),
        if (nzchar(input$sderr)) c("--sderr", input$sderr),
        if (input$fsd != par_def$fsd) c("--fsd", input$fsd),
        if (input$regular_grid == "FALSE") c("--regular_grid", "FALSE"),
        if (nzchar(input$pch) && (input$pch != par_def$pch)) c("--pch", input$pch)
      )
      result <- tryCatch({
        withCallingHandlers({
          do.call(cli, list(args))
        }, warning = function(w) {
          showNotification(paste("Warning:", conditionMessage(w)), type = "warning")
        })
      }, error = function(e) {
        showNotification(paste("Error:", conditionMessage(e)), type = "error")
        c("cli_error", conditionMessage(e))
      })
      args_show=sub(input$meas$datapath, input$meas$name, args, fixed=TRUE)
      args_show=sub(input$sto$datapath, input$sto$name, args_show, fixed=TRUE)
      output$date <- renderText(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
      output$call <- renderText(paste0("R --vanilla -e 'dynafluxr::cli()' --args '", paste0(args_show, collapse="' '"), "'"))
      output$callr <- renderText(paste0("res=dynafluxr::cli(c('", paste0(args_show, collapse="', '"), "'))"))

      if (class(result) == "character" && result[1L] == "cli_error") {
        output$tsv_table <- renderTable(NULL)
        output$pdf_view <- renderUI({
          tags$iframe(style = "width:100%; height:0px;", src = "")
        })
        updateSelectInput(session, "tsv_select", choices = character(0L))
        updateSelectInput(session, "pdf_select", choices = character(0L))
        output$errmes <- renderText(paste0("Error: ", result[2L]))
      } else {
        output$errmes <- renderText(NULL)
        odir=tools::file_path_sans_ext(input$meas$datapath)
        files <- list.files(odir, full.names = TRUE)
        pdf_files <- files[grepl("\\.pdf$", files)]
        tsv_files <- files[grepl("\\.tsv$", files)]
        file.copy(files, temp_dir, overwrite=TRUE)
        if (length(files)) {
          updateSelectInput(session, "pdf_select", choices = basename(pdf_files), selected="ispecie.pdf")
          updateSelectInput(session, "tsv_select", choices = basename(tsv_files), selected="stats.tsv")
          output$tsv_table <- renderTable({
            req(input$tsv_select)
            read.table(file.path(odir, input$tsv_select), header = TRUE, sep = "\t")
          }, digits=4)
          output$pdf_view <- renderUI({
            req(input$pdf_select)
            tags$iframe(style = "width:100%; height:800px;", src = paste0("/pdfs/", input$pdf_select))
          })
          disable("downloadZip")
          enable("pdf_select")
          enable("tsv_select")
        } else {
          enable("downloadZip")
          updateSelectInput(session, "pdf_select", choices = character(0L))
          updateSelectInput(session, "tsv_select", choices = character(0L))
          disable("pdf_select")
          disable("tsv_select")
          output$downloadZip <- downloadHandler(
            filename = function() {
              paste0(tools::file_path_sans_ext(input$meas$name), ".zip")
            },
            content = function(fname) {
              file.copy(file.path(dirname(odir), paste0(basename(odir), ".zip")), fname)
            },
            contentType = "application/zip"
          )
          output$tsv_table <- renderTable(NULL)
          output$pdf_view <- renderUI({
            tags$iframe(style = "width:100%; height:0px;", src = "")
          })
        }
      }
      updateActionButton(inputId = "run",
                 label = "Run",
                 icon = icon("bolt"))

      session$sendCustomMessage(type = "gotop", "gotop: done")
    })
    observeEvent(input$save_params, repeat {
      req(input$save_params)
      if (is.integer(input$save_params))
        break
      save_par <- parseSavePath(volumes, input$save_params)
      write.table(data.frame(Value=inp2val(input)), save_par$datapath, sep = "\t", row.names = TRUE, col.names = TRUE)
      break
    })

    observeEvent(input$load_params, repeat {
      req(input$load_params)
      if (is.integer(input$load_params))
        break
      load_par <- parseFilePaths(volumes, input$load_params)
      params <- read.table(load_par$datapath, sep = "\t", header = TRUE)
      updateNumericInput(session, "knot", value = params["knot", "Value"])
      updateNumericInput(session, "norder", value = params["norder", "Value"])
      updateTextInput(session, "lna", value = params["lna", "Value"])
      updateTextInput(session, "increasing", value = params["increasing", "Value"])
      updateTextInput(session, "decreasing", value = params["decreasing", "Value"])
      updateNumericInput(session, "skip", value = params["skip", "Value"])
      updateCheckboxInput(session, "zip", value = as.logical(params["zip", "Value"]))
      updateCheckboxInput(session, "dls", value = as.logical(params["dls", "Value"]))
      updateCheckboxInput(session, "wsd", value = as.logical(params["wsd", "Value"]))
      updateNumericInput(session, "npi", value = params["npi", "Value"])
      updateTextInput(session, "sf", value = params["sf", "Value"])
      updateTextInput(session, "nosf", value = params["nosf", "Value"])
      updateTextInput(session, "sderr", value = params["sderr", "Value"])
      updateNumericInput(session, "fsd", value = params["fsd", "Value"])
      updateCheckboxInput(session, "regular_grid", value = as.logical(params["regular_grid", "Value"]))
      updateTextInput(session, "pch", value = params["pch", "Value"])
      break
    })
    observeEvent(input$stop, {      
      stopApp()
    })
    session$onSessionEnded(stopApp)
  }
  shinyApp(ui, server, options=list(launch.browser = TRUE))
}
