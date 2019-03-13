On estimation of optimal treatment regimes for maximizing t-year survival probability
R. Jiang, W. Lu, R. Song and M. Davidian
J. R. Statist. Soc. B, 79 (2017), 1165 -- 1185


Description for ACTG175 data: The data contain two files: events.dat and all.dat. The first file contains the survival time information while the second file contains all the covariates information.

For the events.dat file:
Columns are (in order): pidnum, cid, time, trt 

where

pidnum         patient id number
cid	       censoring indicator (1 = failure, 0 = censoring)
time	       time to failure or censoring
trt	       treatment indicator (0 = ZDV only; 1 = ZDV + ddI, 2 = ZDV + Zal, 3 = ddI only)

For the all.dat file:
Columns are (in order):

pidnum age wtkg hemo homo drugs karnof oprior z30 zprior preanti race
gender str2 strat symptom treat offtrt cd40 cd420 cd496 r cd80 cd820

where

pidnum		patient id number
age             age (yrs) at baseline
wtkg            weight (kg) at baseline
hemo  		hemophilia (0=no, 1=yes)
homo            homosexual activity (0=no, 1=yes)
drugs           history of IV drug use (0=no, 1=yes)
karnof          Karnofsky score (on a scale of 0-100)
oprior          Non-ZDV antiretroviral therapy pre-175 (0=no, 1=yes)
z30             ZDV in the 30 days prior to 175 (0=no, 1=yes)
zprior          ZDV prior to 175 (0=no, 1=yes)
preanti         # days pre-175 anti-retroviral therapy
race            race (0=White, 1=non-white)
gender          gender (0=F, 1=M)
str2            antiretroviral history (0=naive, 1=experienced)
strat           antiretroviral history stratification
                 (1='Antiretroviral Naive',2='> 1 but <= 52 weeks of prior
                 antiretroviral therapy',3='> 52 weeks)
symptom         symptomatic indicator (0=asymp, 1=symp)
treat           treatment indicator (0=ZDV only, 1=others)
offtrt          indicator of off-trt before 96+/-5 weeks (0=no,1=yes)
cd40            CD4 at baseline
cd420           CD4 at 20+/-5 weeks
cd496           CD4 at 96+/-5 weeks (=-1 if missing)
r               Missing CD4 at 96+/-5 weeks (0=missing, 1 = observed)
cd80            CD8 at baseline
cd820           CD8 at 20+/-5 weeks


Description of R codes for ACTG175 data analysis: It contains two files: RealData.R and RealData-Functions-v3.R.
You need to run RealData.R to get the results, while RealData-Functions-v3.R contains all the reuired functions. 


Wenbin Lu 
Department of Statistics 
North Carolina State University 
5212 SAS Hall
Box 8203 
Raleigh
NC 27695 
USA 

E-mail: lu@stat.ncsu.edu 
Web: http://www4.stat.ncsu.edu/~lu/
