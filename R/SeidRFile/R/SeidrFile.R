#' @useDynLib SeidRFile, .registration=TRUE
#' @importFrom graphics title
#' @importFrom methods new
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom Rcpp evalCpp

setGeneric("algorithms", function(object) standardGeneric("algorithms"))
setGeneric("nodes", function(object) standardGeneric("nodes"))
setGeneric("supplementary", function(object) standardGeneric("supplementary"))
setGeneric("betweenness", function(object) standardGeneric("betweenness"))
setGeneric("closeness", function(object) standardGeneric("closeness"))
setGeneric("eigenvector", function(object) standardGeneric("eigenvector"))
setGeneric("katz", function(object) standardGeneric("katz"))
setGeneric("strength", function(object) standardGeneric("strength"))
setGeneric("pagerank", function(object) standardGeneric("pagerank"))
setGeneric("reset", function(object) standardGeneric("reset"))
setGeneric("get_chunk", function(object, chunksize = 10000) standardGeneric("get_chunk"))
setGeneric("get_ichunk", function(object, index, chunksize = 10000) standardGeneric("get_ichunk"))

setGeneric("seidr_apply",
           function(X, scores = TRUE, ranks = FALSE, flags = FALSE,
                    supp_str = FALSE, supp_int = FALSE, supp_flt = FALSE,
                    indices = FALSE, FUN, ...) {
             standardGeneric("seidr_apply")
           }
)
setGeneric("seidr_chunked_apply",
           function(X, scores = TRUE, ranks = FALSE, flags = FALSE,
                    supp_str = FALSE, supp_int = FALSE, supp_flt = FALSE,
                    indices = FALSE, chunksize = 10000, FUN, ...) {
             standardGeneric("seidr_chunked_apply")
           }
)
setGeneric("seidr_iapply",
           function(X, scores = TRUE, ranks = FALSE, indices = FALSE,
                    index, FUN, ...) {
             standardGeneric("seidr_iapply")
           }
)
setGeneric("seidr_chunked_iapply",
           function(X, scores = TRUE, ranks = FALSE, indices = FALSE,
                    index, chunksize = 10000, FUN, ...) {
             standardGeneric("seidr_chunked_iapply")
           }
)

#' SeidrFile
#'
#' An R class to access SeidrFiles. It is highly inadvisable for a user to
#' create an instance of this class with \code{new(...)} unless they
#' really know what they are doing. Use at your own risk.
#'
#' @slot sf_ptr External pointer to SeidrFile C++ interface
#' @slot sf_h_ptr External pointer to SeidrFileHeader C++ interface
#' @slot sf_attr List of header attributes
#' @slot path File path
#' @slot i Internal index for source node of current edge
#' @slot j Internal index for target node of current edge
#' @slot k Internal index for the number of edges read
#' @export
setClass("SeidrFile", slots = c("sf_ptr" = "externalptr",
                                "sf_h_ptr" = "externalptr",
                                "sf_attr" = "list",
                                "path" = "character",
                                "i" = "integer",
                                "j" = "integer",
                                "k" = "integer"))

#' SeidrFileFromPath
#'
#' This function creates a new \code{\link{SeidrFile-class}} from a file path. Minimal
#' checks to ensure the file exists and is indeed a gzip compressed
#' file.
#'
#' @param path The file path of the SeidrFile
#' @return A \code{\link{SeidrFile-class}}
#' @seealso \code{\link{SeidrFile-class}}
#' @export
SeidrFileFromPath <- function(path = "") {
  if (! file.exists(path))
  {
    stop("File '", path, "' does not exist")
  }
  i <- file.info(path)
  if (i$isdir)
  {
    stop("File '", path, "' is a directory")
  }

  bp <- file(path, "rb")
  b1 <- readBin(bp, "raw")
  b2 <- readBin(bp, "raw")
  close(bp)

  if (! all(c(b1 == 0x1f, b2 == 0x8b))) {
    stop("File '", path, "' does not look to be gzip compressed")
  }

  p <- normalizePath(path)
  sf_ptr <- SeidrFile__ptr(p)
  sf_h_ptr <- SeidrFileHeader__ptr(sf_ptr)
  sf_attr <- SeidrFileHeader__ptr__getAttr(sf_h_ptr)
  names(sf_attr) <- c("nodes", "edges", "algorithms", "dense", "nsupp",
                      "nsupp_str", "nsupp_int", "nsupp_flt",
                      "pagerank", "closeness", "betweenness",
                      "strength", "eigenvector", "katz")
  sf <- new("SeidrFile", sf_ptr = sf_ptr, sf_h_ptr = sf_h_ptr,
            sf_attr = sf_attr, path = p)
}

#' Print SeidrFile metadata
#'
#' @param object A \code{\link{SeidrFile-class}} object
#' @examples
#' path <- system.file("extdata", "aggregated.sf", package = "SeidRFile")
#' sf <- SeidrFileFromPath(path)
#' show(sf)
#' @export
setMethod("show", signature="SeidrFile", function(object) {

  if (! as.logical(SeidrFile__ptr__isopen(object@sf_ptr))) {
    stop("Connection is invalid")
  }

  cat("SeidrFile             ", object@path, "\n")
  cat("  Nodes:              ", object@sf_attr$nodes, "\n")
  cat("  Edges:              ", object@sf_attr$edges, "\n")
  cat("  Algorithms:         ", object@sf_attr$algorithms, "\n")
  cat("  Storage:            ", ifelse(object@sf_attr$dense, "Dense", "Sparse"), "\n")
  cat("  Supplementary Data: ", object@sf_attr$nsupp, "\n")
  cat("    <int>:            ", object@sf_attr$nsupp_int, "\n")
  cat("    <str>:            ", object@sf_attr$nsupp_str, "\n")
  cat("    <flt>:            ", object@sf_attr$nsupp_flt, "\n")
  cat("  Pagerank:           ", ifelse(object@sf_attr$pagerank, "Yes", "No"), "\n")
  cat("  Closeness:          ", ifelse(object@sf_attr$closeness, "Yes", "No"), "\n")
  cat("  Betweenness:        ", ifelse(object@sf_attr$betweenness, "Yes", "No"), "\n")
  cat("  Strength:           ", ifelse(object@sf_attr$strength, "Yes", "No"), "\n")
  cat("  Eigenvector:        ", ifelse(object@sf_attr$eigenvector, "Yes", "No"), "\n")
  cat("  Katz:               ", ifelse(object@sf_attr$katz, "Yes", "No"), "\n")
})

#' Access node vector in a SeidrFile
#'
#' This function returns the nodes present in the \code{SeidrFile} header.
#' Edge indices correspond to the order in this vector.
#'
#' @param object A \code{\link{SeidrFile-class}} object
#' @return A character vector of node names
#' @examples
#' path <- system.file("extdata", "aggregated.sf", package = "SeidRFile")
#' sf <- SeidrFileFromPath(path)
#' head(nodes(sf))
#' @export
setMethod("nodes", signature="SeidrFile", function(object) {
  return(SeidrFileHeader__ptr__getNodes(object@sf_h_ptr))
})

#' Access supplementary data tags in a SeidrFile
#'
#' This function returns the supplementary data tags present
#' in the \code{SeidrFile} header.
#' Note that these are only the descriptive tags, for the values one
#' still needs to unserialize the edge data.
#'
#' @param object A \code{\link{SeidrFile-class}} object
#' @return A character vector of supplementary data tags.
#' @examples
#' path <- system.file("extdata", "aggregated.sf", package = "SeidRFile")
#' sf <- SeidrFileFromPath(path)
#' supplementary(sf)
#' @export

setMethod("supplementary", signature="SeidrFile", function(object) {
  return(SeidrFileHeader__ptr__getSTags(object@sf_h_ptr))
})

#' Access algorithm names in a SeidrFile
#'
#' This function returns the algorithm names present
#' in the \code{SeidrFile} header.
#'
#' @param object A \code{\link{SeidrFile-class}} object
#' @return A character vector of algorithm names
#' @examples
#' path <- system.file("extdata", "aggregated.sf", package = "SeidRFile")
#' sf <- SeidrFileFromPath(path)
#' algorithms(sf)
#' @export
setMethod("algorithms", signature="SeidrFile", function(object) {
  return(SeidrFileHeader__ptr__getAlgs(object@sf_h_ptr))
})

#' Access centrality data in a SeidrFile
#'
#' These functions access different node centrality metrics (if computed)
#' in a SeidrFile. If they weren't computed \code{numeric(0)} is returned.
#'
#' @param object A \code{\link{SeidrFile-class}} object
#' @return A numeric vector of node centrality values
#' @name centr_dummy
#' @examples
#' path <- system.file("extdata", "aggregated.sf", package = "SeidRFile")
#' sf <- SeidrFileFromPath(path)
#' betweenness(sf)
#'
NULL

#' @rdname centr_dummy
#' @export
setMethod("betweenness", signature="SeidrFile", function(object) {
  return(SeidrFileHeader__ptr__getBetweenness(object@sf_h_ptr))
})

#' @rdname centr_dummy
#' @export
setMethod("closeness", signature="SeidrFile", function(object) {
  return(SeidrFileHeader__ptr__getCloseness(object@sf_h_ptr))
})

#' @rdname centr_dummy
#' @export
setMethod("eigenvector", signature="SeidrFile", function(object) {
  return(SeidrFileHeader__ptr__getEigenvector(object@sf_h_ptr))
})

#' @rdname centr_dummy
#' @export
setMethod("katz", signature="SeidrFile", function(object) {
  return(SeidrFileHeader__ptr__getKatz(object@sf_h_ptr))
})

#' @rdname centr_dummy
#' @export
setMethod("strength", signature="SeidrFile", function(object) {
  return(SeidrFileHeader__ptr__getStrength(object@sf_h_ptr))
})

#' @rdname centr_dummy
#' @export
setMethod("pagerank", signature="SeidrFile", function(object) {
  return(SeidrFileHeader__ptr__getPageRank(object@sf_h_ptr))
})

#' Vectorized apply functions for SeidrFiles
#'
#' These functions apply \code{FUN()} to entire SeidrFiles in a (mostly)
#' vectorized manner. \code{FUN()} is passed a vectorized SeidrFile. The
#' \code{seidr_iapply} and \code{seidr_chunked_iapply} versions apply
#' \code{FUN()} to a single score and/or rank.
#'
#' @param X A \code{\link{SeidrFile-class}} object
#' @param FUN A function to execute
#' @param ... Additional arguments passed to \code{FUN()}
#' @param scores A \code{logical} indicating if scores should be vectorized
#' @param ranks A \code{logical} indicating if ranks should be vectorized
#' @param flags A \code{logical} indicating if flags should be vectorized
#' @param supp_str A \code{logical} indicating if supplementary strings should
#'   be vectorized
#' @param supp_int A \code{logical} indicating if supplementary integers should
#'   be vectorized
#' @param supp_flt A \code{logical} indicating if supplementary floats should be
#'   vectorized
#' @param indices A \code{logical} indicating if edge indices should be
#'   vectorized
#' @param chunksize The size of a single chunk
#' @param index Which score and/or rank to apply \code{FUN()}
#' @return For \code{seidr_apply} and \code{seidr_iapply} the return type is
#'   whatever \code{FUN()} returns.
#'
#'   For \code{seidr_chunked_apply} and \code{seidr_chunked_iapply} a list of
#'   chunks, each containing whatever \code{FUN()} returned for this chunk.
#' @examples
#' # Open a connection to a SeidrFile
#' path <- system.file("extdata", "aggregated.sf", package = "SeidRFile")
#' sf <- SeidrFileFromPath(path)
#'
#' # Get the best rank for each edge
#' best_rank <- function(x) {
#'   return(sapply(x$rank, min, na.rm = TRUE))
#' }
#' best <- seidr_apply(sf, scores = FALSE, ranks = TRUE, FUN = best_rank)
#'
#' # Reset the connection so we can read again
#' sf <- reset(sf)
#'
#' # Create a data.frame from a SeidrFile
#'
#' to_df <- function(x, nodes) {
#'   data.frame("From"=nodes[x$i], "To"=nodes[x$j],
#'              "Weight"=x$score, "Rank"=x$rank)
#' }
#'
#' df <- seidr_iapply(sf, scores = TRUE, ranks = TRUE, indices = TRUE,
#'                    index = 10, FUN = to_df, nodes = nodes(sf))
#'
#' @name apply_dummy
NULL

#' @rdname apply_dummy
#' @export
setMethod("seidr_apply", "SeidrFile",
          function(X, scores = TRUE, ranks = FALSE, flags = FALSE, supp_str = FALSE,
                   supp_int = FALSE, supp_flt = FALSE,
                   indices = FALSE,  FUN, ...) {

            if (! as.logical(SeidrFile__ptr__isopen(X@sf_ptr))) {
              stop("Connection is invalid")
            }

            if(! is.function(FUN)) stop("Not a function: ", deparse(substitute(FUN)))
            pos <- SeidrFileHeader__ptr__Position(X@sf_h_ptr)
            k <- pos[[3]]
            if (k != 0) {
              warning(k, " edges have been read before the call to ",
                      "seidr_apply and will not be included in the ",
                      "output")
            }

            vs <- SeidrFile__ptr__vectorizeSF(X@sf_ptr, X@sf_h_ptr, scores, ranks, flags,
                                              supp_str, supp_int, supp_flt, indices,
                                              chunksize = 0)
            names(vs) <- c("score", "rank", "i", "j", "flag", "supp_str", "supp_int", "supp_flt")
            FUN(vs, ...)

          }
)

#' @rdname apply_dummy
#' @export
setMethod("seidr_chunked_apply", "SeidrFile",
          function(X, scores = TRUE, ranks = FALSE, flags = FALSE, supp_str = FALSE,
                   supp_int = FALSE, supp_flt = FALSE,
                   indices = FALSE, chunksize = 10000, FUN, ...) {

            if (! as.logical(SeidrFile__ptr__isopen(X@sf_ptr))) {
              stop("Connection is invalid")
            }

            if(! is.function(FUN)) stop("Not a function: ", deparse(substitute(FUN)))
            pos <- SeidrFileHeader__ptr__Position(X@sf_h_ptr)
            k <- pos[[3]]
            if (k != 0) {
              warning(k, " edges have been read before the call to ",
                      "seidr_chunked_apply and will not be included in the ",
                      "output")
            }

            nedges <- X@sf_attr$edges
            chunks <- ceiling((nedges - k) / chunksize)
            message("Creating ", chunks, " chunks")
            retl <- vector("list", chunks)
            pb <- txtProgressBar(min = 1, max = chunks, initial = 1, char = "=",
                                 width = NA, title, label, style = 3, file = "")
            for(i in 1:chunks) {
              vs <- SeidrFile__ptr__vectorizeSF(X@sf_ptr, X@sf_h_ptr, scores, ranks, flags,
                                                supp_str, supp_int, supp_flt, indices,
                                                chunksize)
              names(vs) <- c("score", "rank", "i", "j", "flag", "supp_str", "supp_int", "supp_flt")
              retl[[i]] <- FUN(vs, ...)
              setTxtProgressBar(pb, i)
            }
            close(pb)
            return(retl)
          }
)

#' @rdname apply_dummy
#' @export
setMethod("seidr_iapply", "SeidrFile",
          function(X, scores = TRUE, ranks = FALSE, indices = FALSE,
                   index, FUN, ...) {

            if (! as.logical(SeidrFile__ptr__isopen(X@sf_ptr))) {
              stop("Connection is invalid")
            }

            if(! is.function(FUN)) stop("Not a function: ", deparse(substitute(FUN)))
            pos <- SeidrFileHeader__ptr__Position(X@sf_h_ptr)
            k <- pos[[3]]
            if (k != 0) {
              warning(k, " edges have been read before the call to ",
                      "seidr_iapply and will not be included in the ",
                      "output")
            }
            vs <- SeidrFile__ptr__vectorizeSingle(X@sf_ptr, X@sf_h_ptr, scores, ranks, indices,
                                                  index, chunksize = 0)
            names(vs) <- c("score", "rank", "i", "j")
            FUN(vs, ...)
          }
)

#' @rdname apply_dummy
#' @export
setMethod("seidr_chunked_iapply", "SeidrFile",
          function(X, scores = TRUE, ranks = FALSE, indices = FALSE,
                   index, chunksize = 10000, FUN, ...) {

            if (! as.logical(SeidrFile__ptr__isopen(X@sf_ptr))) {
              stop("Connection is invalid")
            }

            if(! is.function(FUN)) stop("Not a function: ", deparse(substitute(FUN)))
            pos <- SeidrFileHeader__ptr__Position(X@sf_h_ptr)
            k <- pos[[3]]
            if (k != 0) {
              warning(k, " edges have been read before the call to ",
                      "seidr_chunked_iapply and will not be included in the ",
                      "output")
            }

            nedges <- X@sf_attr$edges
            chunks <- ceiling((nedges - k) / chunksize)
            message("Creating ", chunks, " chunks")
            retl <- vector("list", chunks)
            pb <- txtProgressBar(min = 1, max = chunks, initial = 1, char = "=",
                                 width = NA, title, label, style = 3, file = "")
            for(i in 1:chunks) {

              vs <- SeidrFile__ptr__vectorizeSingle(X@sf_ptr, X@sf_h_ptr, scores, ranks, indices,
                                                    index, chunksize)
              names(vs) <- c("score", "rank", "i", "j")
              retl[[i]] <- FUN(vs, ...)
              setTxtProgressBar(pb, i)
            }
            close(pb)
            return(retl)
          }
)

#' Close the connection to a SeidrFile
#'
#' @param con A \code{\link{SeidrFile-class}} object
#' @param ... Ignored
#' @export
setMethod("close", signature="SeidrFile", function(con, ...) {
  SeidrFile__ptr__close(con@sf_ptr)
})

#' Regenerate the connection to a SeidrFile
#'
#' This function will reset the SeidrFile for reading.
#'
#' @param object A \code{\link{SeidrFile-class}} object
#' @export
setMethod("reset", "SeidrFile", function(object) {
  warning("Resetting SeidrFile ", deparse(substitute(object)))
  close(object)
  return(SeidrFileFromPath(object@path))
})

#' Retrieve one chunk of edge data from a SeidrFile
#'
#' These functions are mostly meant for development and quick inspection of SeidrFile
#' contents.
#'
#' @param object A \code{\link{SeidrFile-class}} object
#' @param index The index of the score/rank to retrieve
#' @param chunksize The size of the chunk to retrieve
#' @return A list of length \code{chunksize} containing edge data
#' @examples
#' path <- system.file("extdata", "aggregated.sf", package = "SeidRFile")
#' sf <- SeidrFileFromPath(path)
#'
#' get_chunk(sf, chunksize = 10)
#' @name chunk_dummy
NULL

#' @rdname chunk_dummy
#' @export
setMethod("get_chunk", "SeidrFile", function(object, chunksize = 10000) {
  r <- SeidrFile__ptr__vectorizeSF(object@sf_ptr, object@sf_h_ptr, TRUE, TRUE, TRUE, TRUE,
                                   TRUE, TRUE, TRUE, chunksize)
  names(r) <- c("score", "rank", "i", "j", "flag", "supp_str", "supp_int", "supp_flt")
  return(r)
})

#' @rdname chunk_dummy
#' @export
setMethod("get_ichunk", "SeidrFile", function(object, index, chunksize = 10000) {
  r <- SeidrFile__ptr__vectorizeSingle(object@sf_ptr, object@sf_h_ptr, TRUE, TRUE, TRUE,
                                       index, chunksize)
  names(r) <- c("score", "rank", "i", "j")
  return(r)
})



