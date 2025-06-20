###
###

.pkgname <- "BSgenome.Bb31.NCBI.ASM868v2"

.circ_seqs <- c("cp9","cp26","cp32-1","cp32-3","cp32-4","cp32-6","cp32-7","cp32-8","cp32-9")

.onLoad <- function(libname, pkgname)
{
    if (pkgname != .pkgname)
        stop("package name (", pkgname, ") is not ",
             "the expected name (", .pkgname, ")")
    extdata_dirpath <- system.file("extdata", package=pkgname,
                                   lib.loc=libname, mustWork=TRUE)

    ## Make and export BSgenome object.
    bsgenome <- BSgenome(
        organism="Borreliella burgdorferi b31",
        common_name=NA,
        genome="ASM868v2",
        provider="NCBI",
        release_date=NA,
        source_url=NA,
        seqnames=NULL,
        circ_seqs=.circ_seqs,
        mseqnames=NULL,
        seqs_pkgname=pkgname,
        seqs_dirpath=extdata_dirpath
    )

    ns <- asNamespace(pkgname)

    objname <- pkgname
    assign(objname, bsgenome, envir=ns)
    namespaceExport(ns, objname)

    old_objname <- "Bb31"
    assign(old_objname, bsgenome, envir=ns)
    namespaceExport(ns, old_objname)
}

