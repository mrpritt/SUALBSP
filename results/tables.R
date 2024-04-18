library(tidyverse)
library(here)

## this script reproduces the tables from the paper
results=here("results")

## Tables: Lower bounds
##   Here "opt" are the known optimal values of the SBF-2 instances
lb = read_csv(paste0(results,"/lowerbounds.csv"))
lb %>% rowwise() %>% mutate(across(lm1:lm4|lms1dw|lms1sbbf|lmb,~-rd(.x,opt))) %>% group_by(set,alpha) %>% summarize(N=length(set),across(lm1|lms1|lm2|lm3|lm4|lms1sbbf|lms1dw|lmb,~mean(.x,na.rm=T)),lms1sbbf.t=mean(lms1sbbf.t,na.rm=T),lms1dw.t=mean(lms1dw.t)) %>% ungroup() %>% select(-set) %>% as_tibble()

## Tables: Rule-GRASP
##   Here "sub" are the best known values from Scholl (2013), via Scholl (2019, personal communication).
rg = read_csv(paste0(results,"/rule-grasp.csv"))
rg = rg %>% mutate(across(ub10:ub1000,~rd(.x,sub))) %>% group_by(set,alpha) %>% summarize(across(ub10:ub1000|t10:t1000,mean)) %>% as_tibble()
rg %>% bind_rows(rg %>% summarize(across(ub10:last_col(),mean)))

## Tables: Rule sets and Hoffmann
rh = read_csv(paste0(results,"/rules+hoffmann.csv"))

rh %>% separate(inst,c("set","inst"),sep="/") %>% separate(set,c("set","alpha"),sep="-") %>% group_by(set,alpha,alg) %>% summarize(rlb=mean(rd(ub,lb)),time=mean(elapsed),nopt=100*sum(lb==ub)/n()) %>% pivot_wider(names_from=alg,values_from=c(rlb,time,nopt)) %>% as_tibble()

## Tables: CBFS
cbfs = read_csv(paste0(results,"/cbfs.csv"))

cbfs %>% separate(inst,c("set","inst"),sep="/") %>% separate(set,c("set","alpha"),sep="-") %>% group_by(set,alpha,alg) %>% summarize(rlb=mean(rd(ub,lb)),time=mean(elapsed),nopt=100*sum(lb==ub)/n()) %>% pivot_wider(names_from=alg,values_from=c(rlb,time,nopt)) %>% as_tibble()
