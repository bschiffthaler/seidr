library(org.Hs.eg.db)
library(tidyverse)
library(here)

# Get mappings - Entrez ID
## Positive links
p <- AnnotationDbi::select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), columns=c('PATH'))
### Divide by path ID
pmem <- split(p$ENTREZID, p$PATH)
### Filter pathways with kess than 2 members
f <- sapply(pmem, length) > 1
### Combine all possible genes within a PW
r <- t(do.call(cbind, lapply(pmem[f], combn, 2)))
MASS::write.matrix(r, here("goldStandard/hsa_entrez_id_pos.tsv"))

## Negative standard
f <- ! is.na(p$PATH)

### Spread as wide ORF ~ PATH boolean matrix
isec <- as_tibble(p[f, ]) %>% mutate(val = 1) %>% 
  pivot_wider(id_cols = c(ENTREZID, PATH), values_from = val,
              names_from = PATH, values_fill = list(val = 0)) 
### Remove ORF column and cast to actual matrix class
mat <- as.matrix(isec[, -1])
gen <- isec$ENTREZID

### Get all ORF to ORF rowsums
res <- tcrossprod(mat)
dimnames(res) <- list(gen, gen)

### Bind ORFs that have no match in any PATH together as a negative standard
### Could be better vectorized I guess.. fast enough though
neg_std <- do.call(rbind, lapply(1:nrow(res), function(f){
  nomatch <- gen[which(res[f, ] == 0)]
  do.call(rbind, lapply(nomatch, c, gen[f]))
}))

MASS::write.matrix(neg_std, here("goldStandard/hsa_entrez_id_neg.tsv"))

# Get mappings - EnsEMBL
p <- AnnotationDbi::select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), columns=c('PATH','ENSEMBL')) %>% 
  as_tibble() %>% select(ENSEMBL,PATH) %>% unique() %>% filter(!is.na(ENSEMBL)&!is.na(PATH))

## Divide by path ID
pmem <- split(p$ENSEMBL, p$PATH)
## Filter pathways with kess than 2 members
f <- sapply(pmem, length) > 1
## Combine all possible genes within a PW
r <- t(do.call(cbind, lapply(pmem[f], combn, 2)))
MASS::write.matrix(r, here("goldStandard/hsa_ensembl_id_pos.tsv"))

## Negative standard
### Spread as wide ORF ~ PATH boolean matrix
isec <- as_tibble(p) %>% mutate(val = 1) %>% 
  pivot_wider(id_cols = c(ENSEMBL, PATH), values_from = val,
              names_from = PATH, values_fill = list(val = 0)) 

### Remove ORF column and cast to actual matrix class
mat <- as.matrix(isec[, -1])
gen <- isec$ENSEMBL

### Get all ORF to ORF rowsums
res <- tcrossprod(mat)
dimnames(res) <- list(gen, gen)

### Bind ORFs that have no match in any PATH together as a negative standard
### Could be better vectorized I guess.. fast enough though
neg_std <- do.call(rbind, lapply(1:nrow(res), function(f){
  nomatch <- gen[which(res[f, ] == 0)]
  do.call(rbind, lapply(nomatch, c, gen[f]))
}))

MASS::write.matrix(neg_std, here("goldStandard/hsa_ensembl_id_neg.tsv"))

# Get mappins - Symbol
p <- AnnotationDbi::select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), columns=c('PATH','SYMBOL')) %>% 
  as_tibble() %>% select(SYMBOL,PATH) %>% unique() %>% filter(!is.na(SYMBOL)&!is.na(PATH))

## Divide by path ID
pmem <- split(p$SYMBOL, p$PATH)
## Filter pathways with kess than 2 members
f <- sapply(pmem, length) > 1
## Combine all possible genes within a PW
r <- t(do.call(cbind, lapply(pmem[f], combn, 2)))
MASS::write.matrix(r, here("goldStandard/hsa_gene_symbol_pos.tsv"))

## Negative standard
### Spread as wide ORF ~ PATH boolean matrix
isec <- as_tibble(p) %>% mutate(val = 1) %>% 
  pivot_wider(id_cols = c(SYMBOL, PATH), values_from = val,
              names_from = PATH, values_fill = list(val = 0)) 

### Remove ORF column and cast to actual matrix class
mat <- as.matrix(isec[, -1])
gen <- isec$SYMBOL

### Get all ORF to ORF rowsums
res <- tcrossprod(mat)
dimnames(res) <- list(gen, gen)

### Bind ORFs that have no match in any PATH together as a negative standard
### Could be better vectorized I guess.. fast enough though
neg_std <- do.call(rbind, lapply(1:nrow(res), function(f){
  nomatch <- gen[which(res[f, ] == 0)]
  do.call(rbind, lapply(nomatch, c, gen[f]))
}))

MASS::write.matrix(neg_std, here("goldStandard/hsa_gene_symbol_neg.tsv"))

