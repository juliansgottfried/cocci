library(tidyverse)

setwd("~/Desktop/cocci")
variants <- read.table("variantcalls.txt", sep = '\t', header = T)

homozygotes <- as.numeric(sapply(str_split(sapply(str_split(variants$INFO,";"),\(x) x[3]),"="),\(x) x[2]))==0
variants <- variants[homozygotes,]

called <- as.numeric(sapply(str_split(sapply(str_split(variants$INFO,";"),\(x) x[5]),"="),\(x) x[2]))==0
variants <- variants[called,]

polymorphic <- as.numeric(sapply(str_split(sapply(str_split(variants$INFO,";"),\(x) x[4]),"="),\(x) x[2]))!=7
variants <- variants[polymorphic,]

variants <- variants %>% arrange(CHROM, POS)
keep <- rep(F, nrow(variants))
prev_idx <- 1
keep[1] <- T
for (i in 2:nrow(variants)) {
    if (abs(variants$POS[i] - variants$POS[prev_idx]) >= 100) {
        prev_idx <- i
        keep[i] <- T
    }
}
variants <- variants[keep, ]

genotypes <- data.frame(apply(variants %>% select(Sample1:Sample7), 2, function(x) {
    as.numeric(sapply(str_split(sapply(str_split(x, ":"), \(x) x[1]), "/"), \(x) x[1]))
    }))

genotypes <- genotypes %>%
    add_column(.before = "Sample1", Major = variants$REF) %>%
    add_column(.before = "Sample1", Minor = variants$ALT)

for (i in 1:nrow(genotypes)) {
    if (sum(genotypes[i, 3:9], na.rm = T) > 3) {
        genotypes[i, 3:9] = abs(genotypes[i, 3:9] - 1)
        genotypes[i, 1:2] = genotypes[i, 2:1]
    }
}

genotypes <- genotypes %>%
    add_column(.before = "Major", Chromosome = variants$CHROM) %>% 
    add_column(.before = "Major", Position = variants$POS)

genotypes <- na.omit(genotypes)

                