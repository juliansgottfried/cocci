library(tidyverse)

# colors <- c("#492705","#85531c","#c28d46","#d9d1bf")
# colors <- c("#45383f","#755840","#8c5863","#a697ae","#2c8492","#20afef")
# colors <- c("#383434","#587696","#3c9cbb","#fc8f37","#ffb969")
# colors <- c("#45383f","#6e4673","#2c8492","#40bbf5","#f7dcb5")
colors <- c("#45383f","#734650","#2c8492","#40bbf5","#f7dcb5")
bg.color <- colors[1]
bg.color <- "#362c31"
text.color <- colorspace::lighten(colors[1],amount=0.8)

custom_theme <-
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

setwd("~/Desktop/cocci")
variants <- read.table("genomic_data/california_variants.txt", sep = '\t', header = T)
nonSJV <- paste0("Sample",c(1, 5:9))
variants <- variants[,!(colnames(variants) %in% nonSJV)]
colnames(variants)[10:20] <- paste0("Sample",1:11)

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

nalt <- str_split(variants$ALT,",",simplify = T)
keepnonmulti <- sapply(nalt[,2], \(x) x=="")

variants <- variants[keepnonmulti, ]

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
paste0(with(genotypes,ifelse(Sample1,Minor,Major)),collapse="") %>% str_length

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
newdata %>%
    mutate(time = time*scale) %>% 
    ggplot(aes(x = time, y = all)) +
    geom_line(color=colors[4],linewidth = 1) +
    ylim(0, 1) +
    ylab("Ne")+
    theme_classic()+
    custom_theme
newdata %>%
    ggplot(aes(x = time, y = all)) +
    geom_line(color=colors[4],linewidth = 1) +
    ylim(0, 1) +
    ylab("Ne")+
    theme_classic()+
    custom_theme


# some genomic plots
# NW_004504307.1 NW_004504308.1 NW_004504309.1 NW_004504310.1
filtered <- genotypes %>% filter(Chromosome == "NW_004504310.1")
filtered$p <- rowSums(filtered %>% select(Sample1:Sample17)) / 17
filtered <- filtered %>% filter(p <= 0.5)
browser <- read.table("genomic_data/genome_browser.txt", header=T)
# ggplot(filtered, aes(x = p)) +
#     geom_bar() +
#     theme_classic()
chunk_idx <- 2529000
ggplot(filtered, aes(x = Position)) +
    geom_histogram(bins = 50, fill=colors[5], color = colors[2]) +
    xlab("position")+ylab("SNPs")+
    theme_classic()+
    custom_theme
ggplot(browser, aes(x = chromStart)) +
    geom_histogram(bins = 50, fill=colors[2], color = colors[5])+
    xlab("position")+ylab("coding regions")+
    theme_classic()+
    custom_theme
# ggplot(filtered, aes(y = Position, x = as.factor(p))) +
#     geom_violin() +
#     geom_hline(yintercept = 2529000, color = "red") +
#     theme_classic()

min(filtered$Position)

# some other plots
rhohat <- read_csv("outputs/rhohat.csv", col_names = F)
rhopost <- read_csv("outputs/rhohat_posterior.csv", col_names = F)
p0 <- read_csv("outputs/p0.csv", col_names = F)
simrho0 <- read_csv("outputs/sim_est_rho0.csv", col_names = F)
simrho1 <- read_csv("outputs/sim_est_rho1.csv", col_names = F)
r2 <- read_csv("outputs/r2.csv", col_names = F)
baseline <- read_csv("outputs/baseline.csv", col_names = F)

rhos = seq(from=0, to = 19, by = 0.5)

makeaxes <- function(input) {
    input <- input %>% 
        mutate(alpha = seq(from=0, to = 19, by = 0.5)) %>% 
        pivot_longer(cols = -alpha, names_to = "beta")
    input$beta <- rhos[as.numeric(str_split(input$beta, "X",simplify = T)[,2])]
    input
}

rhohat %>% 
    pivot_longer(cols = X1:X2, names_to="X", values_to="Y") %>% 
    ggplot(aes(x=Y, fill=X))+
    labs(fill="parameter", x= "estimate")+
    scale_fill_manual(labels = c("alpha", "beta"), values = c(colors[3], colors[5])) +
    geom_density(color=colors[1], alpha=0.5)+
    theme_classic()+
    custom_theme

p0 <- makeaxes(p0)
ggplot(p0,aes(x=alpha,y=beta,fill=value))+
    geom_tile()+
    scale_fill_gradientn(
        name="log",
        colors=(colors))+
    ggtitle("P(beta = 0)")+
    xlab("alpha")+ylab("beta")+
    theme_classic()+
    coord_equal()+
    custom_theme

thresh <- -235.67767799412619
p0$sig <- p0$value < thresh
p0$sig <- factor(p0$sig, levels=c(T,F))
ggplot(p0,aes(x=alpha,y=beta,fill=sig))+
    geom_tile()+
    scale_fill_discrete(
        name="significant",
        palette=c(colors[5],bg.color),
        labels=c("yes","no"))+
    xlab("alpha")+ylab("beta")+
    theme_classic()+
    coord_equal()+
    custom_theme

rhopost <- makeaxes(rhopost)
#rhopost[38, 3] <- NA
ggplot(rhopost,aes(x=alpha,y=beta,fill=value))+
    geom_tile()+
    scale_fill_gradientn(
        name="log",
        colors=colors, na.value="red")+
    ggtitle("Posterior probability")+
    xlab("alpha")+ylab("beta")+
    theme_classic()+
    coord_equal()+
    custom_theme

simrho0 <- makeaxes(simrho0)
ggplot(simrho0,aes(x=alpha,y=beta,fill=value))+
    geom_tile()+
    scale_fill_gradientn(
        name="estimate",
        colors=colors)+
    xlab("alpha")+ylab("beta")+
    theme_classic()+
    coord_equal()+
    custom_theme

ggplot(simrho0,aes(x=alpha,y=value,color=beta, group=beta))+
    geom_line()+
    scale_color_gradientn(
        name="beta",
        colors=colors[2:5])+
    geom_abline(aes(slope = 1,intercept=0),linetype = 1, color="red")+
    xlab("alpha")+ylab("estimate")+
    theme_classic()+
    coord_equal()+
    xlim(0,19)+ylim(0,19)+
    custom_theme
simrho0 %>% filter(value<1.5) %>% 
    ggplot(aes(x=alpha,y=value,group=beta))+
    geom_line(color="white")+
    geom_abline(aes(slope = 1,intercept=0),linetype = 1, color="red")+
    xlab("alpha")+ylab("estimate")+
    theme_classic()+
    coord_equal()+
    xlim(0,19)+ylim(0,19)+
    custom_theme

simrho1 <- makeaxes(simrho1)
ggplot(simrho1,aes(x=alpha,y=beta,fill=value))+
    geom_tile()+
    scale_fill_gradientn(
        name="estimate",
        colors=colors)+
    xlab("alpha")+ylab("beta")+
    theme_classic()+
    coord_equal()+
    custom_theme
ggplot(simrho1,aes(x=beta,y=value,color=alpha, group=alpha))+
    geom_line()+
    scale_color_gradientn(
        name="alpha",
        colors=colors[2:5])+
    geom_abline(aes(slope = 1,intercept=0),linetype = 1, color="red")+
    xlab("beta")+ylab("estimate")+
    theme_classic()+
    coord_equal()+
    xlim(0,19)+ylim(0,19)+
    custom_theme
simrho1 %>% filter(alpha<4) %>% 
    ggplot(aes(x=beta,y=value, group=alpha))+
    geom_line(color="white")+
    geom_abline(aes(slope = 1,intercept=0),linetype = 1, color="red")+
    xlab("beta")+ylab("estimate")+
    theme_classic()+
    coord_equal()+
    xlim(0,19)+ylim(0,19)+
    custom_theme


r2 <- r2 %>% 
    mutate(dist = seq(from=0, to=1, by=0.01)) %>% 
    pivot_longer(cols = -dist, names_to = "eta")
etas <- rev(seq(from=0.9519,by=0.002,to=0.9999))
r2$eta <- etas[as.numeric(str_split(r2$eta, "X",simplify = T)[,2])]
ggplot(r2,aes(x=dist,y=value,color=1-eta,group=eta))+
    geom_line()+
    scale_color_gradientn(
        name="crossover\nrate",
        colors=((colors[1:4])))+
    xlab("distance")+ylab("r^2")+
    theme_classic()+
    theme(panel.background=element_rect(fill=colors[5],color=colors[5]),
          plot.background=element_rect(fill=colors[5],color=colors[5]),
          legend.background=element_rect(fill=colors[5],color=colors[5]),
          text=element_text(color=colors[1],size=12,family="mono"),
          axis.text=element_text(color=colors[1]),
          axis.text.x=element_text(vjust=0),
          axis.title.y=element_text(vjust=3),
          axis.title.x=element_text(vjust=-1),
          axis.line=element_line(color=colors[1]),
          axis.ticks=element_line(color=colors[1]),
          legend.title=element_text(vjust=4),
          plot.margin=margin(r=15,t=15,l=15,b=15))

ggplot(baseline,aes(x=seq(from=0,by=0.2,to=19.8),y=X2))+
    geom_ribbon(aes(ymin = X1, ymax=X3), fill=colors[3], alpha=0.3)+
    geom_line(color=colors[5])+
    geom_abline(aes(slope = 1,intercept=0),linetype = 2, color="red")+
    xlab("rho")+ylab("estimate")+
    theme_classic()+
    custom_theme


samplingdist <- read_csv("outputs/sampling_distr.csv", col_names = F)
for (i in 1:16) {
    for (j in 1:i) {
        samplingdist[j, i] = samplingdist[i, j]
    }
}
samplingdist <- samplingdist %>% 
    mutate(n1 = 1:16) %>% 
    pivot_longer(cols = -n1, names_to = "n2")
samplingdist$n2 <- as.numeric(str_split(samplingdist$n2, "X",simplify = T)[,2])

ggplot(samplingdist,aes(x=n1,y=n2,fill=log(value/sum(value))))+
    geom_tile()+
    scale_fill_gradientn(
        name="log P",
        colors=colors)+
    xlab("x1")+ylab("x2")+
    theme_classic()+
    coord_equal()+
    ggtitle("x3 = 0")+
    theme(
          text=element_text(size=12,family="mono"),
          axis.text.x=element_text(vjust=0),
          axis.title.y=element_text(vjust=3),
          axis.title.x=element_text(vjust=-1),
          legend.title=element_text(vjust=4),
          plot.margin=margin(r=15,t=15,l=15,b=15))


rhoprob <- read_csv("outputs/rhoprob.csv", col_names = F)
rhoprob <- rhoprob %>% rename(probability = "X1") %>% mutate(rho = seq(from=0, to = 19.9, by = 0.1))
ggplot(rhoprob,aes(x=rho,y=probability))+
    geom_line()+
    theme_classic()+
    theme(
        text=element_text(size=12,family="mono"),
        axis.text.x=element_text(vjust=0),
        axis.title.y=element_text(vjust=3),
        axis.title.x=element_text(vjust=-1),
        legend.title=element_text(vjust=4),
        plot.margin=margin(r=15,t=15,l=15,b=15))



data.frame(x = 1:100, y = 1) %>% 
    ggplot(aes(x,y))+
    geom_line()+
    theme_classic()+
    labs(x="time", y="rate")+
    ylim(0,2)+
    theme(
        text=element_text(size=12,family="mono"),
        axis.text.x=element_text(vjust=0),
        axis.title.y=element_text(vjust=3),
        axis.title.x=element_text(vjust=-1),
        legend.title=element_text(vjust=4),
        plot.margin=margin(r=15,t=15,l=15,b=15))
data.frame(x = 1:100, y = 0.5*(2+sin(exp(25:124/30)))) %>% 
    ggplot(aes(x,y))+
    geom_line()+
    theme_classic()+
    labs(x="time", y="rate")+
    ylim(0,2)+
    theme(
        text=element_text(size=12,family="mono"),
        axis.text.x=element_text(vjust=0),
        axis.title.y=element_text(vjust=3),
        axis.title.x=element_text(vjust=-1),
        legend.title=element_text(vjust=4),
        plot.margin=margin(r=15,t=15,l=15,b=15))



alphaq <- quantile(rhohat$X1,c(0.055,0.5,0.945))
betaq <- quantile(rhohat$X2,c(0.055,0.5,0.945))
betaribbon <- newdata %>% 
    select(time, all) %>% 
    mutate(time = scale*time) %>% 
    mutate(upper = betaq[3]*all, lower = betaq[1]*all, rho = betaq[2]*all, component = "beta")
alpharibbon <- newdata %>% 
    select(time, all) %>% 
    mutate(time = scale*time) %>% 
    mutate(upper = alphaq[3], lower = alphaq[1], rho = alphaq[2],  component = "alpha")
ribbons <- rbind(alpharibbon, betaribbon)
ribbons %>% ggplot(aes(time, rho, color=component))+
    geom_line()+
    geom_ribbon(aes(ymin=lower,ymax=upper, fill=component), alpha=0.2, color=F)+
    scale_color_manual(labels = c("alpha", "beta * M(t)"), values = c(colors[4], colors[5])) +
    scale_fill_manual(guide="none",values = c(colors[4], colors[5])) +
    ylim(0,19)+
    theme_classic()+
    custom_theme

loglik <- read_csv("outputs/loglik.csv", col_names = F)
data.frame(L=loglik$X1, rho = rhos) %>% 
    ggplot(aes(x=rho,y=L))+
    geom_line()+
    ylab("log lik")+
    geom_vline(aes(xintercept = rhos[which.max(loglik$X1)]),color="red",linetype=2)+
    theme_classic()+
    theme(
        text=element_text(size=12,family="mono"),
        axis.text.x=element_text(vjust=0),
        axis.title.y=element_text(vjust=3),
        axis.title.x=element_text(vjust=-1),
        legend.title=element_text(vjust=4),
        plot.margin=margin(r=15,t=15,l=15,b=15))
