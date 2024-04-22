sum_for_T <- function(X1, X2)
{
  i = 1
  
  sum = 0
  
  n1 = length(rownames(X1))
  n2 = length(rownames(X2))
  
  X1_transpose = t(X1)
  
  while(i < n1){
    
    j = 1
    
    while(j < n2){
      if(identical(X1, X2))
      {
        if(i != j)
        {
          sum = sum + X1_transpose[1,i]*X2[j,1]
        }
      }
      else
      {
        sum = sum + X1_transpose[1,i]*X2[j,1]
      }
      
      j = j+1
    }
    
    i = i+1
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