#' Translate plain text file with reactions to stoichiometric matrix
#' 
#' @param fn String, file name
#' @return Matrix of size n_metab x n_reac
#' @importFrom utils read.delim
#' @export

txt2sto=function(fn) {
    te=read.delim(fn, header=FALSE, row.names=1L, comment.char="#")
    sto=data.frame()
    for (re in row.names(te)) {
        #if (re == "vSUC")
        #    browser()
        #lr=trimws(strsplit(te[re,], "->", fixed=TRUE)[[1L]])
        lr=trimws(strsplit(paste(te[re,], sep="", collapse=NULL), "->", fixed=TRUE)[[1L]])
        lr[[1]]=sub(" *<*-*$", "", lr[[1]])
        for (ip in 1:2) {
            if (ip > length(lr))
                next
            part=lr[ip]
            mets=trimws(strsplit(part, "+", fixed=TRUE)[[1L]])
            for (pr in mets) {
                terms=trimws(strsplit(pr, "*", fixed=TRUE)[[1L]])
                s=(ip-1.5)*2
                if (length(terms) == 1L) {
                    coef=s
                    met=terms
                } else {
                    coef=s*as.double((terms[1]))
                    met=terms[2]
                }
                if (is.null(sto[met, re]) || is.na(sto[met, re])) {
                    sto[met, re]=coef
                    sto[is.na(sto)]=0.
                } else {
                    sto[met, re]=sto[met, re] + coef
                }
            }
        }
    }
    as.matrix(sto)
}
