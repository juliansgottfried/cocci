library(tidyverse)

setwd("~/Desktop/cocci")
variants <- read.table("genomic_data/california_variants.txt", sep = '\t', header = T)

n <- length(grep("Sample",colnames(variants)))

# variants <- variants %>% filter(CHROM != "NW_027094050.1")

homozygotes <- as.numeric(sapply(str_split(sapply(str_split(variants$INFO,";"),\(x) x[3]),"="),\(x) x[2]))==0
variants <- variants[homozygotes,]

called <- as.numeric(sapply(str_split(sapply(str_split(variants$INFO,";"),\(x) x[5]),"="),\(x) x[2]))==0
variants <- variants[called,]

polymorphic <- as.numeric(sapply(str_split(sapply(str_split(variants$INFO,";"),\(x) x[4]),"="),\(x) x[2]))!=n
variants <- variants[polymorphic,]

# variants <- variants %>% arrange(CHROM, POS)
# keep <- rep(F, nrow(variants))
# prev_idx <- 1
# keep[1] <- T
# for (i in 2:nrow(variants)) {
#     if (abs(variants$POS[i] - variants$POS[prev_idx]) >= 100) {
#         prev_idx <- i
#         keep[i] <- T
#     }
# }
# variants <- variants[keep, ]

genotypes <- data.frame(apply(variants %>% select(starts_with("Sample")), 2, function(x) {
    as.numeric(sapply(str_split(sapply(str_split(x, ":"), \(x) x[1]), "/"), \(x) x[1]))
    }))

genotypes <- genotypes %>%
    add_column(.before = "Sample1", Major = variants$REF) %>%
    add_column(.before = "Sample1", Minor = variants$ALT)

for (i in 1:nrow(genotypes)) {
    if (sum(genotypes[i, 3:(n+2)], na.rm = T) > n%/%2) {
        genotypes[i, 3:(n+2)] = abs(genotypes[i, 3:(n+2)] - 1)
        genotypes[i, 1:2] = genotypes[i, 2:1]
    }
}

genotypes <- genotypes %>%
    add_column(.before = "Major", Chromosome = variants$CHROM) %>% 
    add_column(.before = "Major", Position = variants$POS)

genotypes <- na.omit(genotypes)

write_csv(genotypes,"genomic_data/california_data.csv")

# get whole sample sequences
paste0(with(genotypes,ifelse(Sample7,Minor,Major)),collapse="")


# get constant sites
chroms <- unique(genotypes$Chromosome)
getconst <- function(i) {
    snps <- c("A"=0,"C"=0,"G"=0,"T"=0)
    name <- paste0("genomic_data/venezuela_chrom",i,".txt")
    chrom <- paste0(read.table(name)$V1,collapse="")
    tmpsnps <- genotypes$Major[genotypes$Chromosome == chroms[i]] %>% table
    snps[names(snps) %in% names(tmpsnps)] = tmpsnps
    str_count(chrom, names(snps)) - snps
}
apply(sapply(1:8, getconst), 1, sum)

