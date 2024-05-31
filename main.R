library(ALL)
library(geneset)
library(genefilter)
library(tidyverse)

import::from(sum_for_T_stat.R, sum_for_T)
import::from(trace_for_Q_stat.R, trace_for_Q)

data(ALL)

# кількість генів
# genes amount
print("genes:")
print(length(row.names(exprs(ALL))))

# кількість пацієнтів
# patients amount
print("patients:")
print(length(colnames(exprs(ALL))))

# from Gentleman et al. "Bioinformatics and Computational Biology Solutions Using R and Bioconducto
subset <- intersect(grep("^B", as.character(ALL$BT)),
                    + which(ALL$mol %in% c("BCR/ABL", "NEG")))
ALL <- ALL[, subset]

# from Gentleman et al. "Bioinformatics and Computational Biology Solutions Using R and Bioconductor
# прибрати некеспресовані та низькомінливі гени
# remove not expressed and low variability genes
f1 <- pOverA(0.25, log2(100))
f2 <- function(x) (IQR(x) > 0.5)
ff <- filterfun(f1, f2)
selected <- genefilter(ALL, ff)
sum(selected)
ALL <- ALL[selected, ]

# кількість генів після прибирання
# genes amount after removal
print("genes:")
print(length(row.names(exprs(ALL))))

# кількість пацієнтів після відбору
# patients amount after selection
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

# прибирає негенні колонки
# removing non-gene columns
X1 = X1[,1:(ncol(X1)-21)]
X2 = X2[,1:(ncol(X2)-21)]

X1_means = lapply(X=X1, FUN=mean)

X2_means = lapply(X=X2, FUN=mean)

sample_genes_amount = 10

set.seed(1)

samplegenes = sample(colnames(X1), sample_genes_amount)

colors = c('black', 'red')

plotsymbol = 19

x = seq(1, length(row.names(exprs(ALL))))

plot(x, X1_means, pch = plotsymbol, col = colors[1], cex = 1.5, xlab = 'genes', ylab = 'expression level', type = 'l')

title(main = 'Gene expression level means BCR/ABL')

plot(x, X2_means, pch = plotsymbol, col = colors[2], cex = 1.5, type = 'l', xlab = 'genes', ylab = 'expression level')

title(main = 'Gene expression level means NEG')


sign_level = 0.05

start <- Sys.time()

sum_X1 <- sum_for_T(X1, X1)

sum_X2 <- sum_for_T(X2, X2)

sum_X1X2 <- sum_for_T(X1, X2)

T_n <- sum_X1 + sum_X2 + sum_X1X2

trace_X1 <- trace_for_Q(X1, X1)

trace_X2 <- trace_for_Q(X2, X2)

trace_X1X2 <- trace_for_Q(X1, X2)


sigma_n <- (2/(n_1*(n_1-1)))*trace_X1 + (2/(n_2*(n_2-1)))*trace_X2 + (4/(n_1*n_2))*trace_X1X2

Q_n <- T_n / sqrt(sigma_n)

print("T_n = ")

print(T_n)

print("Q_n = ")

print(Q_n)

Q_norm = qnorm(1-sign_level/2)

print("Q_norm")

print(Q_norm)

print("Execution time")

print(Sys.time()-start)

if(Q_n > Q_norm)
{
  print("Means are not equal")
} else
{
  print("Means are equal")
}
