trace_for_Q <- function(X1, X2)
{
  # потрібно прописати обчислення
  
  trace = 1
  
  n1 = length(rownames(X1))
  n2 = length(rownames(X2))
  
  if(identical(X1, X2)) {
    trace = trace*(2/(n1*(n1-1))) 
  }
  else{
    trace = trace*(4/(n1*n2)) 
  }
  
  return(trace)
  
}