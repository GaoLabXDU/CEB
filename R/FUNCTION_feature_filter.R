# FUNCTION_feature_filter

filter_params <- function(dataset, fraction) {
  n.cells <- dim(dataset)[2]
  
  min.cells <- ceiling(fraction*n.cells)
  max.cells <- ceiling(fraction*n.cells)
  min.reads <- 2
  
  return(list(min.cells = min.cells, max.cells = max.cells,
              min.reads = min.reads))
}

feature_filter <- function(data, fraction) {
  cat("Feature filtering...\n")
  filter.params <- filter_params(data, fraction)
  min.cells <- filter.params$min.cells
  max.cells <- filter.params$max.cells
  min.reads <- filter.params$min.reads
  d <- data[rowSums(data > min.reads) >= min.cells &
              rowSums(data > 0) <= dim(data)[2] - max.cells, ]
  d <- unique(d)
  if(dim(d)[1] == 0) {
    cat("All features were removed after the feature filter! Stopping now...")
    return()
  }
  return(d)
}