#' @export
encode_mrfi <- function(vec_included){
  paste0(as.numeric(vec_included), collapse = "")
}

#' @export
decode_mrfi <- function(coded_char){
  as.logical(as.integer(strsplit(coded_char, split = "")[[1]]))
}
