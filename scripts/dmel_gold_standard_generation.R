library(org.Dm.eg.db)
library(tidyverse)
library(here)
# Get mappings
p <- AnnotationDbi::select(org.Dm.eg.db, keys=keys(org.Dm.eg.db), columns=c('PATH', 'FLYBASE'))
# Divide by path ID
pmem <- split(p$FLYBASE, p$PATH)
# Filter pathways with kess than 2 members
f <- sapply(pmem, length) > 1
# Combine all possible genes within a PW
r <- t(do.call(cbind, lapply(pmem[f], combn, 2)))
MASS::write.matrix(r, here("goldStandard/dmel_pos.tsv"))


#### Negative standard
f <- ! is.na(p$PATH)
# Spread as wide ORF ~ PATH boolean matrix
isec <- as_tibble(p[f, c(3, 2)]) %>% mutate(val = 1) %>% 
  pivot_wider(id_cols = c(FLYBASE, PATH), values_from = val,
              names_from = PATH, values_fill = list(val = 0)) 
# Remove ORF column and cast to actual matrix class
mat <- as.matrix(isec[, -1])
gen <- isec$FLYBASE
# Get all ORF to ORF rowsums
res <- tcrossprod(mat)
dimnames(res) <- list(gen, gen)
# Bind ORFs that have no match in any PATH together as a negative standard
# Could be better vectorized I guess.. fast enough though
neg_std <- do.call(rbind, lapply(1:nrow(res), function(f){
  nomatch <- gen[which(res[f, ] == 0)]
  do.call(rbind, lapply(nomatch, c, gen[f]))
}))

MASS::write.matrix(neg_std, here("goldStandard/dmel_neg.tsv"))

