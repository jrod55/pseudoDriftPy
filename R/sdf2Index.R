#' @title sdf2index
#' @description Ensure data is in valid sdf format then index the sdf file. sdf files can be downloaded from https://mona.fiehnlab.ucdavis.edu/downloads. \cr
#' This function makes use of \link[ChemmineR]{} functions to ensure validity and formatting is correct for pseudoDrift
#' @name sdf2index
#' @param sdf_input Input file in .sdf format.
#' @param out_name Name of the .xls file to be written out.
#' @return Two files written to the current working directory: 1) an index .xls file with the \code{outname} prefix, 2) an sdf file with the \code{valid} prefix.
#' @importFrom ChemmineR read.SDFset
#' @importFrom ChemmineR validSDF
#' @importFrom ChemmineR write.SDF
#' @importFrom ChemmineR datablock2ma
#' @importFrom ChemmineR datablock
#' @importFrom ChemmineR sdfStream
#' @export
#' @examples
#' sdf2index(sdf_input = system.file("extdata", "test.sdf", package = "pseudoDrift"), out_name = "Index.xls")
sdf2index <- function(sdf_input, out_name){
  # inputs check
  if (!grepl("\\.sdf$", sdf_input)) {
    stop(paste0("Check your input file name: ", sdf_input), " is not a .sdf file")
  }
  if (!grepl("\\.xls$", out_name)) {
    stop(paste0("Check your output file name: ", out_name), " needs to end in .xls")
  }
  sdfset <- suppressWarnings(ChemmineR::read.SDFset(paste0(sdf_input)))
  valid <- ChemmineR::validSDF(sdfset)
  sdfset <- sdfset[valid]
  sdf_input_name = basename(sdf_input)
  ChemmineR::write.SDF(sdfset, file=paste0("valid-",sdf_input_name))
  desc <- function(sdfset) {
    cbind(
      ChemmineR::datablock2ma(datablocklist = ChemmineR::datablock(x = sdfset))[,c("NAME","ID")]
    )
  }
  ChemmineR::sdfStream(input = paste0("valid-",sdf_input_name),
                       output = paste0(out_name),
                       fct=desc, Nlines=1e3, silent = TRUE)
}
