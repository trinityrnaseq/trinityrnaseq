source("../misc_rnaseq_funcs.R")
source("../heatmap.3.R")

data = mtcars;
data = log2(data*10+1)
scale(data,center=T,scale=T); # Z-scores


## Use principal components to partition the 'genes'
pc = princomp(data, cor=T)

pc4 = pc$scores[,1:4]
pc4_cols = matrix_to_color_assignments(pc4)
colnames(pc4_cols) = paste("PC", 1:4)

#heatmap.3(data, RowSideColors=pc4_cols)


## Highlight samples 
Pon = matrix(grepl("c", colnames(data)), nrow=1)
print(Pon)
PonMat = apply(Pon, 1:2, function(x) as.numeric(x))
PonMat = PonMat+1

PonMatCol = matrix_to_color_assignments(PonMat, col=c('white','black'))
rownames(PonMatCol) = "Pon"

heatmap.3(data, RowSideColors=pc4_cols, ColSideColors=PonMatCol)
