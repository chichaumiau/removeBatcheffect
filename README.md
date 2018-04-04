# removeBatcheffect
python scripts to remove batch effect

This function is exactly the same as removeBatchEffect function in [limma](https://bioconductor.org/packages/release/bioc/html/limma.html)

limma(pheno, exprs, covariate_formula, design_formula='1', rcond=1e-8):

## parameters:
pheno: the metadata as a table

exprs: the expression profile, the rows are the genes and the columns are the cells.

covariate_formula: is an expression of the variance you want to keep. (columns in the metadata table).

design_formula: is an expression of the variance you want to remove. i.e. the batch. 

## example:
regressed=limma(pheno, exprs, "~ age + cancer -1", "C(batch)")
