# FUNCTION_data_trans


datatype_trans <- function(inputTags, datatype){
  if (datatype == "count"){  # transformed to CPM
      input_transformed <- t(t(inputTags)/colSums(inputTags))*1000000
  } else if (datatype == "CPM" || datatype == "FPKM" || datatype == "RPKM" || datatype == "TPM") {
      input_transformed <- inputTags
  }
}

log_trans <- function(dataset){
      dataset <- log2(dataset + 1)
}