library(ALL)
library(geneset)

import::from(sum_for_T_stat.R, sum_for_T)
import::from(trace_for_Q_stat.R, trace_for_Q)

# if (!require("BiocManager", quietly = TRUE))
  
#  install.packages("BiocManager")

# BiocManager::install("ALLMLL")

data(ALL)

print(ALL)

# кількість генів
print("genes:")
print(length(row.names(exprs(ALL))))

# кількість пацієнтів
print("patients:")
print(length(colnames(exprs(ALL))))

ALL1 <- data.frame(ALL)

# print(colnames(ALL1))

# print(length(colnames(ALL1)))

ALL_BCRABL = ALL1[ALL$mol.biol == 'BCR/ABL', ]

n_1 = length(ALL_BCRABL$X1000_at);

print("BCR/ABL patients:")
print(n_1)

# print(ALL_BCRABL$X1000_at)

ALL_NEG = ALL1[ALL$mol.biol == 'NEG', ]

n_2 = length(ALL_NEG$X1000_at);

print("NEG patients:")
print(n_2)

# print(ALL_NEG$X1000_at)

i = 1

while (i<length(ALL_BCRABL)){
  
  X1 = ALL_BCRABL[i]
  X2 = ALL_NEG[i]
  
  sign_level = 95
  
  sum_x1 <- sum_for_T(X1, X1)

  sum_x2 <- sum_for_T(X2, X2)
  
  sum_X1X2 <- sum_for_T(X1, X2)
  
  T_n = sum_X1 + sum_X2 + sum_X1X2
  
  
  
  trace_X1 <- trace_for_Q(X1, X1)
  
  trace_X2 <- trace_for_Q(X2, X2)
  
  trace_X1X2 <- trace_for_Q(X1, X2)
  
  sigma_n = trace_X1 + trace_X2 + trace_X1X2
  
  Q_n = T_n / sigma_n
   
  print("T_n = ")
  
  print(T_n)
  
  print("Q_n = ")
  
  print(Q_n)
  
  # якщо Q_n > xi_a, де xi_a верхній sign_level квантиль розподілу N(0,1) 
  # (d-вимірного? як знайти d?), то відкидається H_0, середні не рівні
  
  print(" ") 
  
  i = i+1
}

# bp = getGO(data_dir = tempdir(), ont="bp")
# print(summary(bp))
# 
# print(colnames(bp$geneset))
# 
# print(bp$geneset$gene)
# 
# print(length(unique(unlist(bp$geneset$bp, use.names = FALSE))))
# 
# mf = getGO(data_dir = tempdir(), ont="mf")
# print(summary(mf))
# 
# cc = getGO(data_dir = tempdir(), ont="cc")
# print(summary(cc))
# 
# print(cc$geneset)