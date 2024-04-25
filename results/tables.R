library(tidyverse)
library(here)
library(glue)

## this script reproduces the tables from the paper
results=here("results")

normal_names = function(df) {
    df %>% separate(inst,c("set","inst"),sep="/") %>% separate(set,c("set","alpha"),sep="-") %>% mutate(alpha=as.numeric(alpha))
}

######################################################################
## Tables
######################################################################

## Tables: Lower bounds
##   Here "opt" are the known optimal values of the SBF-2 instances
lb = read_csv(paste0(results,"/lowerbounds.csv"))
lb %>% rowwise() %>% mutate(across(lm1:lm4|lms1dw|lms1sbbf|lmb,~-rd(.x,opt))) %>% group_by(set,alpha) %>% summarize(N=length(set),across(lm1|lms1|lm2|lm3|lm4|lms1sbbf|lms1dw|lmb,~mean(.x,na.rm=T)),lms1sbbf.t=mean(lms1sbbf.t,na.rm=T),lms1dw.t=mean(lms1dw.t)) %>% ungroup() %>% select(-set) %>% as_tibble()

## Tables: Rule-GRASP
##   Here "sub" are the best known values from Scholl (2013), via Scholl (2019, personal communication).
rg = read_csv(paste0(results,"/rule-grasp.csv"))
rt = rg %>% mutate(across(ub10:ub1000,~rd(.x,sub))) %>% group_by(set,alpha) %>% summarize(across(ub10:ub1000|t10:t1000,mean)) %>% as_tibble()
rt %>% bind_rows(rg %>% summarize(across(ub10:last_col(),mean)))

## Tables: Rule sets and Hoffmann
rh = read_csv(paste0(results,"/rules+hoffmann.csv"))

rh %>% normal_names() %>% group_by(set,alpha,alg) %>% summarize(rlb=mean(rd(ub,lb)),time=mean(elapsed),nopt=100*sum(lb==ub)/n()) %>% pivot_wider(names_from=alg,values_from=c(rlb,time,nopt)) %>% as_tibble()

## Tables: CBFS
cbfs = read_csv(paste0(results,"/cbfs.csv"))

cbfs %>% normal_names() %>% group_by(set,alpha,alg) %>% summarize(rlb=mean(rd(ub,lb)),time=mean(elapsed),nopt=100*sum(lb==ub)/n()) %>% pivot_wider(names_from=alg,values_from=c(rlb,time,nopt)) %>% as_tibble()

######################################################################
## Additional information from the paper
######################################################################
sch = read_table(paste0(results,"/Scholl,etal (2013).dat")) %>% normal_names()
lfs = read_csv(paste0(results,"/l-f-star.csv"))
rgs = rg %>% full_join(sch)

ub = rg %>% inner_join(lfs) %>% rowwise() %>% mutate(lmb=max(c_across(lm1:lm4),na.rm=T)) %>% select(-c(lm1:lm4)) %>% mutate(nopt = sum(ub1000 == lmb)) %>% group_by(set,alpha) %>% summarize(opt=100*sum(nopt)/n())

## Lower bounds
lm4.agree = 100 * sum(with(lb, lm4==opt)) / nrow(lb)
dw.notagree = sum(with(lb, lms1dw != opt))
cat(glue("LM4 is clearly the best of the simple bounds and agrees with the optimal value on {format(lm4.agree,digits=1)}% of the instances.\n\n"))
cat(glue("LMDW produces the best bound, and on all except {dw.notagree} larger instances equals the optimal value.\n\n"))

## Upper bounds
or1 = mean(ub$opt)
or2 = with(ub %>% filter(set == "SBF1" & alpha >= 0.5), max(opt))
or3 = with(ub %>% filter(!(set == "SBF1" & alpha >= 0.5)), min(opt))

cat(glue("When comparing the upper bounds to L∗f , we find average optimality rates of about {or1}%. Notably, optimality rates for instances SBF1 with α ≥ 0.5 never exceed {or2}% while all other lie above {or3}%.\n\n"))

ub1000.agree = with(rgs, sum(ub1000==ub))
ub1000.more = with(rgs, sum(ub1000>ub))
ub1000.less = with(rgs, sum(ub1000<ub))
ub1000.time = mean(rt$t1000)
cat(glue("We finally remark that FBR1000 agrees on {ub1000.agree} instances with the best found values of Scholl et al. (2013) over all tests. In {ub1000.more} instances it needs an extra station and in {ub1000.less} instances one station less. Note that this is achieved by a single run with an average running time of less than {format(round(ub1000.time,digits=1))} seconds.\n\n"))
