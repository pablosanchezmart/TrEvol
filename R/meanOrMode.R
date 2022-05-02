#' Compute mean if variable is numeric or mode if variable is character or factor.
#'
#' @param x numeric, character or factor vector.
#'
#' @return mean or mode value
#' @export
#'
#' @importFrom rlang .data
meanOrMode <- function(x){
  if(is.numeric(x) | is.double(x) | is.integer(x)){
    meanX <- mean(x, na.rm = T)
    return(meanX)
  }
  if(is.character(x) | is.factor(x)){
    tble <- as.data.frame(table(x))

    modeX <- tble %>% dplyr::filter(.data$Freq == max(.data$Freq)) %>% dplyr::pull(.data$x)

    return(modeX)
  }
}
