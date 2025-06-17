# Function to mask column names
mask_column <- function(x, old_name, new_name){
  if(old_name != new_name){

    # Temporary name if mask already exists
    if(new_name %in% colnames(x)) {
      colnames(x)[colnames(x) == new_name] <- paste0(new_name, "_renamed")
      colnames(x)[colnames(x) == old_name] <- new_name

    # Remove temporary name
    } else {
      colnames(x)[colnames(x) == old_name] <- new_name
      colnames(x)[colnames(x) == paste0(old_name, "_renamed")] <- old_name
    }
  }

  x
}
