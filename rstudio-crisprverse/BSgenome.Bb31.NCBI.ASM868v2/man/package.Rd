\name{BSgenome.Bb31.NCBI.ASM868v2}
\docType{package}

\alias{BSgenome.Bb31.NCBI.ASM868v2-package}
\alias{BSgenome.Bb31.NCBI.ASM868v2}
\alias{Bb31}

\title{Full genomic sequences for Borreliella burgdorferi b31 (NCBI assembly ASM868v2)}

\description{
  Full genomic sequences for Borreliella burgdorferi b31 as provided by NCBI (assembly ASM868v2, assembly accession GCF_000008685.2). The sequences are stored in DNAString objects.
}

\author{Michael Foster <nope@no.co>}

\seealso{
  \itemize{
    \item \link[BSgenome]{BSgenome} objects in the \pkg{BSgenome}
          software package.

    \item The \code{\link[GenomeInfoDb]{seqinfo}} getter and
          \link[GenomeInfoDb]{Seqinfo} objects in the \pkg{GenomeInfoDb}
          package.

    \item The \code{\link[GenomeInfoDb]{seqlevelsStyle}} getter and
          setter in the \pkg{GenomeInfoDb} package.

    \item \link[Biostrings]{DNAString} objects in the \pkg{Biostrings}
          package.

    \item The \code{\link[BSgenome]{available.genomes}} function
          in the \pkg{BSgenome} software package.

    \item The \link[BSgenomeForge]{forgeBSgenomeDataPkgFromNCBI}
          function in the \pkg{BSgenomeForge} package that was used to
          create this package.
  }
}

\examples{
BSgenome.Bb31.NCBI.ASM868v2
bsg <- BSgenome.Bb31.NCBI.ASM868v2
head(seqlengths(bsg))
seqinfo(bsg)

## ---------------------------------------------------------------------
## Genome-wide motif searching
## ---------------------------------------------------------------------

## See the GenomeSearching vignette in the BSgenome software
## package for some examples of genome-wide motif searching using
## Biostrings and the BSgenome data packages:
if (interactive())
    vignette("GenomeSearching", package="BSgenome")
}

\keyword{package}
\keyword{data}
