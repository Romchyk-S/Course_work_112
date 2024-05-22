library(ALL)
library(geneset)
library(genefilter)
library(tidyverse)

import::from(sum_for_T_stat.R, sum_for_T)
import::from(trace_for_Q_stat.R, trace_for_Q)

data(ALL)

# кількість генів
print("genes:")
print(length(row.names(exprs(ALL))))

# кількість пацієнтів
print("patients:")
print(length(colnames(exprs(ALL))))

# from Gentleman et al. "Bioinformatics and Computational Biology Solutions Using R and Bioconducto
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

# print(summary(ALL))

# кількість генів
print("genes:")
print(length(row.names(exprs(ALL))))

# кількість пацієнтів
print("patients:")
print(length(colnames(exprs(ALL))))

ALL1 <- data.frame(ALL)

ALL_BCRABL = ALL1[ALL$mol.biol == 'BCR/ABL', ]

n_1 = length(ALL_BCRABL$X1005_at);

print("BCR/ABL patients:")
print(n_1)

ALL_NEG = ALL1[ALL$mol.biol == 'NEG', ]

n_2 = length(ALL_NEG$X1005_at);

print("NEG patients:")
print(n_2)

X1 = ALL_BCRABL
X2 = ALL_NEG

# прибирає негенні колонки.
X1 = X1[,1:(ncol(X1)-21)]
X2 = X2[,1:(ncol(X2)-21)]

print(summary(t(X1)))

print("")

print(summary(t(X2)))

# set.seed(1)

samplegenes = sample(colnames(X1), 6)

colors = c('black', 'red')

plotsymbol = 19

for(gene1 in samplegenes)
{
  for(gene2 in samplegenes)
  {
    if(!identical(gene1, gene2))
    {
      plot(ALL_BCRABL[,gene1], ALL_BCRABL[,gene2], pch = plotsymbol, col = colors[1], cex = 1.5, xlab = gene1, ylab = gene2)
      points(ALL_NEG[,gene1], ALL_NEG[,gene2], pch = plotsymbol, cex = 1.5, col = colors[2])
      
      legend(x = "topleft", title = "Leukemia classes", legend = c("BCR/ABL", "NEG"), col = c(colors[1], colors[2]), cex = 0.75, pch = plotsymbol)
    }
  }
}

# 
# sign_level = 0.05
# 
# start <- Sys.time()
# 
# sum_X1 <- sum_for_T(X1, X1)
# 
# sum_X2 <- sum_for_T(X2, X2)
# 
# sum_X1X2 <- sum_for_T(X1, X2)
# 
# T_n <- sum_X1 + sum_X2 + sum_X1X2
# 
# trace_X1 <- trace_for_Q(X1, X1)
# 
# trace_X2 <- trace_for_Q(X2, X2)
# 
# trace_X1X2 <- trace_for_Q(X1, X2)
# 
# 
# sigma_n <- (2/(n_1*(n_1-1)))*trace_X1 + (2/(n_2*(n_2-1)))*trace_X2 + (4/(n_1*n_2))*trace_X1X2
# 
# Q_n <- T_n / sqrt(sigma_n)
# 
# print("T_n = ")
# 
# print(T_n)
# 
# print("Q_n = ")
# 
# print(Q_n)
# 
# Q_norm = qnorm(1-sign_level/2)
# 
# print("Q_norm")
# 
# print(Q_norm)
# 
# print("Execution time")
# 
# print(Sys.time()-start)
# 
# if(Q_n > Q_norm)
# {
#   print("Means are not equal")
# } else
# {
#   print("Means are equal")
# }

# якщо Q_n > xi_a, де xi_a верхній sign_level квантиль розподілу N(0,1)
#, то відкидається H_0, середні не рівні
