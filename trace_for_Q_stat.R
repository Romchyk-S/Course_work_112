library(zoo)
library(dplyr)
library(psych)

trace_for_Q <- function(X1, X2)
{
  trace = 0
  
  n1 = length(rownames(X1))
  n2 = length(rownames(X2))
  
  X1 <- as.data.frame(X1)
  
  X2 <- as.data.frame(X2)
  
  for(i in rownames(X1))
  {
    for(j in rownames(X2))
    {
      X1_vector = as.numeric(c(X1[i,]))
      
      X2_vector = as.numeric(c(X2[j,]))
      
      if(identical(X1,X2))
      {
        X1_copy = data.matrix((X1[!(row.names(X1) %in% c(i,j)),]))
      
        X2_copy = data.matrix((X2[!(row.names(X2) %in% c(i,j)),]))
      }
      else
      {
        X1_copy = data.matrix((X1[!(row.names(X1) %in% c(i)),]))
        
        X2_copy = data.matrix((X2[!(row.names(X2) %in% c(j)),]))
      }
      
      X1_copy = na.spline(X1_copy)
      
      X2_copy = na.spline(X2_copy)
      
      sample_mean_1 = mean(X1_copy)

      sample_mean_2 = mean(X2_copy)

      if(!(identical(X1, X2)) | i != j)
      {
        X1_vector = data.matrix(X1_vector)

        X2_vector = data.matrix(X2_vector)

        result_matrix = (X1_vector-sample_mean_1)%*%t(X1_vector)%*%(X2_vector-sample_mean_2)%*%t(X2_vector)
        
        trace = trace + tr(result_matrix)
      }
    }
  }
  
  if(identical(X1, X2))
  { 
    trace = trace / (n1*(n1-1))
    
  }
  else
  {
    trace = trace/(n1*n2)
  }
  
  return(trace)
  
}