sum_for_T <- function(X1, X2)
{
  j = 1
  
  sum = 0
  
  n1 = length(rownames(X1))
  n2 = length(rownames(X2))
  
  while(j < n1){
    
    k = 1
    
    while(k < n2){
      if(identical(X1, X2))
      {
        if(j != k)
        {
          sum = sum + X1[j,1]*X2[k,1]
        }
      }
      else
      {
        sum = sum + X1[j,1]*X2[k,1]
      }
      
      k = k+1
    }
    
    j = j+1
  }
  
  if(identical(X1, X2))
  {
    
    
    sum = sum / (n1*(n1-1))
  }
  else
  {
    sum = -2*sum/(n1*n2)
  }
  
  return(sum)
  
}