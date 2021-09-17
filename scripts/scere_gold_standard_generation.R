library(org.Sc.sgd.db)
library(tidyverse)
library(here)
# Get mappings
p <- AnnotationDbi::select(org.Sc.sgd.db, keys=keys(org.Sc.sgd.db), columns=c('PATH'))
# Divide by path ID
pmem <- split(p$ORF, p$PATH)
# Filter pathways with kess than 2 members
f <- sapply(pmem, length) > 1
# Combine all possible genes within a PW
r <- t(do.call(cbind, lapply(pmem[f], combn, 2)))
MASS::write.matrix(r, here("goldStandard/scere_pos.tsv"))

#### Negative standard
f <- ! is.na(p$PATH)
# Spread as wide ORF ~ PATH boolean matrix
isec <- as_tibble(p[f,c(1,3)]) %>% mutate(val = 1) %>% 
  pivot_wider(id_cols = c(ORF, PATH), values_from = val,
              names_from = PATH, values_fill = list(val = 0)) 
# Remove ORF column and cast to actual matrix class
mat <- as.matrix(isec[, -1])
gen <- isec$ORF
# Get all ORF to ORF rowsums
res <- tcrossprod(mat)
dimnames(res) <- list(gen, gen)
# Bind ORFs that have no match in any PATH together as a negative standard
# Could be better vectorized I guess.. fast enough though
neg_std <- do.call(rbind, lapply(1:nrow(res), function(f){
  nomatch <- gen[which(res[f, ] == 0)]
  do.call(rbind, lapply(nomatch, c, gen[f]))
}))

MASS::write.matrix(neg_std, here("goldStandard/scere_neg.tsv"))


# Cell cycle
# Galactose
# Starch and Sucrose
idx <- match(c('04111', '00052', '00500'), names(pmem))
pmemx <- pmem[-idx]
f <- sapply(pmemx, length) > 1
# Combine all possible genes within a PW
r <- t(do.call(cbind, lapply(pmemx[f], combn, 2)))
MASS::write.matrix(r,
                   here("goldStandard/scere_pos_lo.tsv"),
                   sep = '\t')
pmemy <- pmem[idx]
f <- sapply(pmemy, length) > 1
# Combine all possible genes within a PW
r <- lapply(pmemy[f], combn, 2)
r <- lapply(names(f), function(x) {rbind(r[[x]], rep(x, ncol(r[[x]])))} )
r <- t(do.call(cbind, r))
MASS::write.matrix(r,
                   here("goldStandard/scere_pos_lo_test.tsv"),
                   sep = '\t')
