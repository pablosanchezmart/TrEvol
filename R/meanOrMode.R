#' Return mean if variable is numeric or mode if variable is character or factor.
#'
#' @param x (*numeric*, *character* or *factor*). Vector from which the mean or mode will be computed for.
#'
#' @return (*numeric* for mean or *character* for mode). Mean or mode value.
#' @export
#'
#' @importFrom rlang .data
#'
#' #' @examples
#' \dontrun{
#' # Calculate mean
#'
#' numeric.data <- c(1, 2, 2, 3)
#' meanOrMode(numeric.data)
#'
#' # Calculate mean
#'
#' categorical.data <- c("1", "2", "2", "3")
#' meanOrMode(categorical.data)
#' }
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
