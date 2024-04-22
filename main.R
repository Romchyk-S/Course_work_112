library(ALL)
library(geneset)
library(genefilter)

import::from(sum_for_T_stat.R, sum_for_T)
import::from(trace_for_Q_stat.R, trace_for_Q)

# if (!require("BiocManager", quietly = TRUE))
  
#  install.packages("BiocManager")

# BiocManager::install("ALLMLL")

# BiocManager::install("genefilter")

data(ALL)

subset <- intersect(grep("^B", as.character(ALL$BT)),
                    + which(ALL$mol %in% c("BCR/ABL", "NEG")))
ALL <- ALL[, subset]

# from Gentleman et al. "Bioinformatics and Computational Biology Solutions Using R and Bioconductor"
# remove not expressed and low variability genes
f1 <- pOverA(0.25, log2(100))
f2 <- function(x) (IQR(x) > 0.5)
ff <- filterfun(f1, f2)
selected <- genefilter(ALL, ff)
sum(selected)
ALL <- ALL[selected, ]

# кількість генів
print("genes:")
print(length(row.names(exprs(ALL))))

# print(row.names(exprs(ALL)))

# кількість пацієнтів
print("patients:")
print(length(colnames(exprs(ALL))))

ALL1 <- data.frame(ALL)

# print(colnames(ALL1))

# print(length(colnames(ALL1)))

ALL_BCRABL = ALL1[ALL$mol.biol == 'BCR/ABL', ]

n_1 = length(ALL_BCRABL$X1005_at);

print("BCR/ABL patients:")
print(n_1)

ALL_NEG = ALL1[ALL$mol.biol == 'NEG', ]

n_2 = length(ALL_NEG$X1005_at);

print("NEG patients:")
print(n_2)


i = 1

sum_X1 = 0
sum_X2 = 0
sum_X1X2 = 0

while (i<length(ALL_BCRABL)){

  X1 = ALL_BCRABL[i]
  X2 = ALL_NEG[i]

  sign_level = 0.05
  
  # Bonferroni correction

  sum_x1 <- sum_for_T(X1, X1)

  sum_x2 <- sum_for_T(X2, X2)

  sum_X1X2 <- sum_for_T(X1, X2)

  T_n = sum_X1 + sum_X2 + sum_X1X2



  trace_X1 <- trace_for_Q(X1, X1)

  trace_X2 <- trace_for_Q(X2, X2)

  trace_X1X2 <- trace_for_Q(X1, X2)
  

  sigma_n = (2/(n_1*(n_1-1)))*trace_X1 + (2/(n_2*(n_2-1)))*trace_X2 + (4/(n_1*n_2))*trace_X1X2
  
  print(sigma_n)
  
  print(sqrt(sigma_n))

  Q_n = T_n / sqrt(sigma_n)

  print("T_n = ")

  print(T_n)

  print("Q_n = ")

  print(Q_n)

  # якщо Q_n > xi_a, де xi_a верхній sign_level квантиль розподілу N(0,1)
  #, то відкидається H_0, середні не рівні

  print(" ")

  i = i+1
  
  # break
}

bp = getGO(data_dir = tempdir(), ont="bp")
print(summary(bp))

# print(row.names(bp$geneset_name))

# print(bp$geneset$gene)

# print(bp$geneset$gene)

# print(length(unique(unlist(bp$geneset$bp, use.names = FALSE))))
# 
# mf = getGO(data_dir = tempdir(), ont="mf")
# print(summary(mf))
# 
# print(length(unique(unlist(mf$geneset$mf, use.names = FALSE))))
# 
# cc = getGO(data_dir = tempdir(), ont="cc")
# print(summary(cc))
# 
# print(length(unique(unlist(cc$geneset$cc, use.names = FALSE))))

tykin <- unique(lookup("GO:0004713", "hgu95av2",
                       + "GO2ALLPROBES"))
print(tykin)

# print(cc$geneset)