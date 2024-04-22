trace_for_Q <- function(X1, X2)
{
  # потрібно прописати обчислення
  
  trace = 0
  
  n1 = length(rownames(X1))
  n2 = length(rownames(X2))
  
  if(identical(X1, X2)) {
    
    j = 1
    
    k = 1
    
    while(j<n1)
    {
      while(k<n2)
      {
        
        if(j != k)
        {
          X1_transpose = t(X1)

          X1copy <- X1[-c(1,j),]
        
          X1copy <- X1copy[-c(k)]
          
          sample_mean = mean(X1copy)
          
          trace = trace + (X1[j,1]-sample_mean)*X1_transpose[1,j]*(X1[k,1]-sample_mean)*X1_transpose[1,k]
        }

        k = k+1
        
      }
      
      j = j+1
    }
    
    print(trace)
    
    trace = trace/(n1*(n1-1))
  }
    
  else{
    
    l = 1
    
    k = 1
    
    while(l<n1)
    {
      while(k<n2)
      {
        
        X1_transpose = t(X1)
        
        X2_transpose = t(X2)

        X1copy <- X1[-c(1,l),]
        
        X2copy <- X2[-c(1,k),]
        
        sample_mean_1 = mean(X1copy)
        
        sample_mean_2 = mean(X2copy)
        
        trace = trace + (X1[l,1]-sample_mean_1)*X1_transpose[1,l]*(X2[k,1]-sample_mean_2)*X2_transpose[1,k]
    
        k = k+1
        
      }
      
      l = l+1
    }
    
    trace = trace/(n1*n2)
  }
  
  return(trace)
  
}