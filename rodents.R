# install.packages("portalr")
# use_default_data_path("~/Desktop/cocci/")
library(portalr)
library(tidyverse)

data_tables <- load_rodent_data("repo")
download_observations(".")
rodent_data <- abundance(".", time = "date", type = "granivores") %>% 
    mutate(year=year(as_date(censusdate))) %>% 
    pivot_longer(-c(censusdate,year),names_to="species",values_to="pop") %>% 
    group_by(year,species) %>% 
    summarize(pop=mean(pop))

rodent_data %>% 
    ggplot(aes(x=year,y=pop,color=species,group=species))+
    geom_line()+
    theme_classic()

rodent_data %>% 
    summarize(pop=sum(pop)) %>% 
    ggplot(aes(x=year,y=pop))+
    geom_line()+
    theme_classic()
