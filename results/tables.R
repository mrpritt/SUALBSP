library(tidyverse)
library(here)
library(glue)
library(EnvStats)

## this script reproduces the tables from the paper
results=here("results")

normal_names = function(df) {
    df %>% separate(inst,c("set","inst"),sep="/") %>% separate(set,c("set","alpha"),sep="-") %>% mutate(alpha=as.numeric(alpha))
}
rd = function(v,o) { 100*(v-o)/o }

######################################################################
## Tables
######################################################################

## Tables: Lower bounds
##   Here "opt" are the known optimal values of the SBF-2 instances
shift = 1
lb = read_csv(paste0(results,"/lowerbounds.csv"))
lb %>% rowwise() %>% mutate(across(lm1:lm4|lms1dw|lms1sbbf|lmb,~-rd(.x,opt))) %>% group_by(set,alpha) %>% summarize(N=length(set),across(lm1|lms1|lm2|lm3|lm4|lms1sbbf|lms1dw|lmb,~mean(.x,na.rm=T)),lms1sbbf.t=geoMean(lms1sbbf.t+shift,na.rm=T)-shift,lms1dw.t=geoMean(lms1dw.t+shift)-shift) %>% ungroup() %>% select(-set) %>% as_tibble()

## Tables: Rule sets and Hoffmann
rh = read_csv(paste0(results,"/rules+hoffmann.csv"))
rh %>% normal_names() %>% group_by(set,alpha,alg) %>% summarize(rlb=mean(rd(ub,lb)),time=mean(elapsed),nopt=100*sum(lb==ub)/n()) %>% pivot_wider(names_from=alg,values_from=c(rlb,time,nopt)) %>% as_tibble()

## Tables: CBFS
cbfs = read_csv(paste0(results,"/cbfs.csv"))
cbfs %>% normal_names() %>% group_by(set,alpha,alg) %>% summarize(rlb=mean(rd(ub,lb)),time=mean(elapsed),nopt=100*sum(lb==ub)/n()) %>% pivot_wider(names_from=alg,values_from=c(rlb,time,nopt)) %>% as_tibble()
