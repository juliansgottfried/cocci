library(tidyverse)
x <- seq(from=0,to=10,by=0.01)
y <- 1-exp(-0.5*x)
data.frame(x=x,y=y) %>% 
    ggplot(aes(x,y)) +
    geom_line() +
    labs(x="Time",y="",title="Probability of death 1-e⁻") +
    scale_x_continuous(breaks = 0:10, labels = c(0,rep("",9),"T"))+
    theme_classic()+
    theme(text=element_text(size=18,family="mono"))
