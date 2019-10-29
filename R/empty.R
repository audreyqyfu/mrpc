#Check if a data frame is empty

empty <- function (df)
{
  (is.null(df) || nrow(df) == 0 || ncol(df) == 0)
}
