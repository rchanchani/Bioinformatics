source("http://tinyurl.com/rescale-R")
rescale
rescale( c(1,10,1:1000))
rescale(1:10)

## write a both_na() function
# start with a simple setup where we know what the answer should be
x <- c(1,2,NA,3,NA)
y <- c(NA,3,NA,3,4)
 
is.na(x)
is.na(y)

both_na3 <- function(x, y) {
  ## Print some info on where NA's are as well as the number of them 
  if(length(x) != length(y)) {
    stop("Input x and y should be vectors of the same length", call.=FALSE)
  }
  na.in.both <- ( is.na(x) & is.na(y) )
  na.number  <- sum(na.in.both)
  na.which   <- which(na.in.both)
  
  message("Found ", na.number, " NA's at position(s):", 
          paste(na.which, collapse=", ") ) 
  
  return( list(number=na.number, which=na.which) )
}

