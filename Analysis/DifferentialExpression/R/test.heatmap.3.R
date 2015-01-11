
# EXAMPLE USAGE

source("heatmap.3.R")
 
# example of colsidecolors rowsidecolors (single column, single row)
mat <- matrix(1:100, byrow=T, nrow=10)
column_annotation <- sample(c("red", "blue", "green"), 10, replace=T)
column_annotation <- matrix(column_annotation, nrow=1)
rownames(column_annotation) <- c("Variable X")
 
row_annotation <- sample(c("red", "blue", "green"), 10, replace=T)
row_annotation <- matrix(row_annotation, ncol=1)
colnames(row_annotation) <- c("Variable Y")
 
heatmap.3(mat, RowSideColors=row_annotation, ColSideColors=column_annotation)

cat("\n","Press Enter to Continue","\n")
readLines(file("stdin"),1)

#######################
# two annotations each 
#######################

mat <- matrix(1:100, byrow=T, nrow=10)
column_annotation <- matrix(sample(c("red", "blue", "green"), 20, replace=T), nrow=2)
rownames(column_annotation) <- c("Variable X1", "Variable X2")
 
row_annotation <- matrix(sample(c("red", "blue", "green"), 20, replace=T), ncol=2)
colnames(row_annotation) <- c("Variable Y1", "Variable Y2")
 
heatmap.3(mat, RowSideColors=row_annotation, ColSideColors=column_annotation)
 
cat("\n","Press Enter to Continue","\n")
readLines(file("stdin"),1)

########################
# four annotations each
########################

mat <- matrix(1:100, byrow=T, nrow=10)
column_annotation <- matrix(sample(c("red", "blue", "green"), 40, replace=T), nrow=4)
rownames(column_annotation) <- c("Var X1", "Var X2", "Var X3", "Var X4")

row_annotation <- matrix(sample(c("red", "blue", "green"), 40, replace=T), ncol=4)
colnames(row_annotation) <- c("Var Y1", "Var Y2", "Var Y3", "Var Y4")

heatmap.3(mat, RowSideColors=row_annotation, ColSideColors=column_annotation)

cat("\n","Press Enter to Continue","\n")
readLines(file("stdin"),1)

