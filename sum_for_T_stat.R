library(zoo)
library(dplyr)

sum_for_T <- function(X1, X2)
{
  sum = 0
  
  n1 = length(rownames(X1))
  n2 = length(rownames(X2))

  for(i in rownames(X1))
  {

    for(j in rownames(X2))
    {
      X1_vector = as.numeric(c(X1[i,]))
      
      X2_vector = as.numeric(c(X2[j,]))
      
      X1_vector = na.spline(X1_vector)
      
      X2_vector = na.spline(X2_vector)
      
      if(!(identical(X1, X2)) | i != j)
      {
        sum = sum + t(X1_vector)%*%X2_vector
      }
    }
  }
  
  print(sum)
  
  if(identical(X1, X2))
  { 
    sum = sum / (n1*(n1-1))
    
  }
  else
  {
    sum = -2*sum/(n1*n2)
  }

  print(sum)
  
  print(" ")

  
}