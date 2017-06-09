
expect_vectors_equal <- function(a,b){
  if (sum(a %in% b) == length(a) & length(setdiff(a,b)) == 0){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

expect_has_all_and_only_these_members <- function(l, v){
  return(expect_vectors_equal(names(l), v))
}
