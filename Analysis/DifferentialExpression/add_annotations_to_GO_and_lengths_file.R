

# contributed by: Satyajeet Khare  satyajeetkhare@gmail.com
# May 3, 2018

# If you modified your gene count matrix or transcript count matrix to add functional annotations using ...
# Trinity/Analysis/DifferentialExpression/rename_matrix_feature_identifiers.pl script, you will get an error
# if you run Trinity/Analysis/DifferentialExpression/analyze_diff_expr.pl script for Gene Ontology analysis.
# The error will be such as "Error in gene_lengths[features_with_GO, ] : subscript out of bounds".
# You will get this error because of descrepancy between Trinity gene ids in "go_annotations" file and ...
# annotated gene IDs in "gene.matrix" file in edger or deseq out folder.
# To correct this error, modify the go_annotation file and gene_length file using the script below.

setwd("Your_working_directory/")

# Create data object for Trinity.gene.lengths
Trinity_gene_lengths <- read.csv("Trinity.gene_lengths.txt", sep = "\t")

# Create data object for go_annotations
go_annotations <- read.csv("go_annotations.txt", sep = NULL, header = FALSE)

# Create data object for annot_feature_map
annot_feature_map <- read.csv("annot_feature_map.txt", sep = "\t", header = FALSE)

# Modify Trinity_gene_lengths data object by adding a third column with Annotations
Trinity_gene_lengths_mod <- left_join(Trinity_gene_lengths, annot_feature_map, by = c("X.gene_id" = "V1"))

# Modify third column by replacing NA with Trinity gene IDs
setDT(Trinity_gene_lengths_mod)[is.na(V2), V2 := X.gene_id]

# Relace the "X.gene_id" values with values in the third column
Trinity_gene_lengths_mod[, "X.gene_id"] <- Trinity_gene_lengths_mod$V2

# Delete the third column. Its no longer required.
Trinity_gene_lengths_mod$V2 <- NULL

# Write a modified gene length file. This file will be used for analysis of differential expression
write.table(Trinity_gene_lengths_mod, file = "Trinity_gene_lengths_mod.txt", quote = FALSE, row.names = FALSE, sep = "\t")

# Modify go_annotations by adding a third column with Annotations
go_annotations_mod <- left_join(go_annotations, annot_feature_map, by = c("V1" = "V1"))

# Relace the "V1" values with values in the third column (V2.y)
go_annotations_mod[, "V1"] <- go_annotations_mod$V2.y

# Delete the third column. Its no longer required.
go_annotations_mod$V2 <- NULL

# Write a modified go annotation file. This file will be used for analysis of differential expression
write.table(go_annotations_mod, file = "go_annotations_mod.txt", quote = FALSE, row.names = FALSE, sep = "\t")

# Use "Trinity_gene_lengths_mod.txt" in place of "Trinity.gene_lengths.txt" 
# and "go_annotations_mod.txt" in place of "go_annotations.txt"
# in "Trinity/Analysis/DifferentialExpression/analyze_diff_expr.pl"script.
