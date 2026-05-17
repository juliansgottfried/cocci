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

# estimate Ne
mu_bp <- 10^(-10)

L <- 28.9 * 10^6
mu <- mu_bp*L
Ne <- (nrow(genotypes)/(2*sum(1/(1:(n-1)))))/mu
Ne <- 10^7

# estimate gen time in years
tmrca <- 2.85*10^(-5)
bounds <- c(1500, 5500)
gs <- bounds * mu_bp / tmrca

# get the scaling!
scales <- Ne*gs

# write data
chrom <- genotypes$Chromosome %>% table %>% which.max %>% names
filtered <- genotypes %>% filter(Chromosome == chrom)
filtered %>% 
    select(Sample1:Sample17, Position) %>% 
    mutate(Position = Position - min(Position) + 1) %>% 
    write_csv(col_names = F, "sampling_data/alleles.csv")

# smooth rodent data
rodents <- read_csv("rodent_data/rodents.csv")
sp <- colnames(rodents)[2:ncol(rodents)]
rodents$time <- rodents$time / scales[2]
nsp <- ncol(rodents) - 1
dt <- 0.01
mintime <- ceiling(min(rodents$time) / dt) * dt
newdata = data.frame(time = seq(from = mintime, 
                                to = max(rodents$time), 
                                by = 0.01))
for (i in 1:nsp) {
    fit <- loess(pull(rodents[, i + 1]) ~ time, rodents, span = 0.25)
    newdata[sp[i]] <- predict(fit, newdata)
}
newdata$all <- apply(newdata[, 2:(nsp + 1)], 1, sum)
newdata[, 2:(nsp + 2)] <- apply(newdata[, 2:(nsp + 2)], 2, \(x) x / max(x))
adder <- cbind(data.frame(time = seq(from = dt, by = dt, to = mintime - dt)),
      newdata[1, 2:ncol(newdata)])
newdata <- rbind(adder, newdata)
newdata %>% 
    pivot_longer(cols = -time, names_to = "species", values_to = "Ne") %>% 
    ggplot(aes(x = time, y = Ne)) +
    geom_line() +
    ylim(0, 1) +
    xlim(0, 3) +
    facet_wrap(~species, ncol = 1) +
    theme_classic()
write_csv(newdata %>% select(time, stephensi_chartreuse), col_names = F, "rodent_data/covariate.csv")
