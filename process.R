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

# estimate scale
mu_bp_gen <- 10^(-10)
mu_bp_yr <- 1.08 * 10^(-9)
# L <- 28.9 * 10^6
# mu <- mu_bp_gen*L
# Ne <- (nrow(genotypes)/(2*sum(1/(1:(n-1)))))/mu
Ne <- 10^7
g <- mu_bp_gen / mu_bp_yr
# bounds <- c(1500, 5500)
# tmrca <- 2.85*10^(-5)
# g <- bounds * mu_bp_gen / tmrca
scale <- Ne*g

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
rodents$time <- rodents$time / scale
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
    facet_wrap(~species, ncol = 1) +
    theme_classic()
write_csv(newdata %>% select(time, all), col_names = F, "rodent_data/covariate.csv")


# some genomic plots
# NW_004504307.1 NW_004504308.1 NW_004504309.1 NW_004504310.1
filtered <- genotypes %>% filter(Chromosome == "NW_004504310.1")
filtered$p <- rowSums(filtered %>% select(Sample1:Sample17)) / 17
filtered <- filtered %>% filter(p <= 0.5)
ggplot(filtered, aes(x = p)) +
    geom_bar() +
    theme_classic()
ggplot(filtered, aes(x = Position)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = 2529000, color = "red") +
    theme_classic()
ggplot(filtered, aes(y = Position, x = as.factor(p))) +
    geom_violin() +
    geom_hline(yintercept = 2529000, color = "red") +
    theme_classic()

# some other plots
rhohat <- read_csv("outputs/rhohat.csv", col_names = F)
rhopost <- read_csv("outputs/rhohat_posterior.csv", col_names = F)
p0 <- read_csv("outputs/p0.csv", col_names = F)

coquina <- c("#1e0b0c","#8c618f","#b0aba5","#dec885","white")
bg.color <- coquina[1]
text.color <- colorspace::lighten(coquina[1],amount=0.8)

rhohat %>% 
    pivot_longer(cols = X1:X2, names_to="X", values_to="Y") %>% 
    ggplot(aes(x=Y, fill=X))+
    labs(fill="parameter", x= "estimate")+
    scale_fill_manual(labels = c("alpha", "beta"), values = c(coquina[2], coquina[4])) +
    geom_density(color=coquina[1], alpha=0.7)+
    theme_classic()+
    theme(panel.background=element_rect(fill=bg.color,color=bg.color),
          plot.background=element_rect(fill=bg.color,color=bg.color),
          legend.background=element_rect(fill=bg.color,color=bg.color),
          text=element_text(color=text.color,size=12,family="mono"),
          axis.text=element_text(color=text.color),
          axis.text.x=element_text(vjust=0),
          axis.title.y=element_text(vjust=3),
          axis.title.x=element_text(vjust=-1),
          axis.line=element_line(color=text.color),
          axis.ticks=element_line(color=text.color),
          legend.title=element_text(vjust=4),
          plot.margin=margin(r=15,t=15,l=15,b=15))

rhos = seq(from=0, to = 19, by = 0.5)
p0 <- p0 %>% 
    mutate(alpha = seq(from=0, to = 19, by = 0.5)) %>% 
    pivot_longer(cols = -alpha, names_to = "beta")
p0$beta <- rhos[as.numeric(str_split(p0$beta, "X",simplify = T)[,2])]

ggplot(p0,aes(x=alpha,y=beta,fill=value))+
    geom_tile()+
    scale_fill_gradientn(
        name="log P(beta = 0)",
        colors=coquina)+
    xlab("alpha")+ylab("beta")+
    theme_classic()+
    coord_equal()+
    theme(panel.background=element_rect(fill=bg.color,color=bg.color),
          plot.background=element_rect(fill=bg.color,color=bg.color),
          legend.background=element_rect(fill=bg.color,color=bg.color),
          text=element_text(color=text.color,size=12,family="mono"),
          axis.text=element_text(color=text.color),
          axis.text.x=element_text(vjust=0),
          axis.title.y=element_text(vjust=3),
          axis.title.x=element_text(vjust=-1),
          axis.line=element_line(color=text.color),
          axis.ticks=element_line(color=text.color),
          legend.title=element_text(vjust=4),
          plot.margin=margin(r=15,t=15,l=15,b=15))

rhopost <- rhopost %>% 
    mutate(alpha = seq(from=0, to = 19, by = 0.5)) %>% 
    pivot_longer(cols = -alpha, names_to = "beta")
rhopost$beta <- rhos[as.numeric(str_split(rhopost$beta, "X",simplify = T)[,2])]

ggplot(rhopost,aes(x=alpha,y=beta,fill=value))+
    geom_tile()+
    scale_fill_gradientn(
        name="log P",
        colors=coquina)+
    xlab("alpha")+ylab("beta")+
    theme_classic()+
    coord_equal()+
    theme(panel.background=element_rect(fill=bg.color,color=bg.color),
          plot.background=element_rect(fill=bg.color,color=bg.color),
          legend.background=element_rect(fill=bg.color,color=bg.color),
          text=element_text(color=text.color,size=12,family="mono"),
          axis.text=element_text(color=text.color),
          axis.text.x=element_text(vjust=0),
          axis.title.y=element_text(vjust=3),
          axis.title.x=element_text(vjust=-1),
          axis.line=element_line(color=text.color),
          axis.ticks=element_line(color=text.color),
          legend.title=element_text(vjust=4),
          plot.margin=margin(r=15,t=15,l=15,b=15))
