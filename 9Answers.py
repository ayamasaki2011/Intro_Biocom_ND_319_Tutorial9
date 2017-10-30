import pandas
file=pandas.read_csv("ponzr1.csv",header=0,sep=",")

mut1=file.loc[file.mutation.isin(['WT', 'M124K']),:]
mut2=file.loc[file.mutation.isin(['WT', 'V456D']),:]
mut3=file.loc[file.mutation.isin(['WT', 'I213N']),:]

#mutation 1 subset into a new dataframe
mut12=pandas.DataFrame({'y':mut1.ponzr1Counts, 'x':0})
mut12.loc[mut1.mutation=='M124K', 'x']=1

#mutation 2 subset into a new dataframe
mut22=pandas.DataFrame({'y':mut2.ponzr1Counts, 'x':0})
mut22.loc[mut1.mutation=='V456D', 'x']=1

#mutation 3 subset into a new dataframe
mut32=pandas.DataFrame({'y':mut3.ponzr1Counts, 'x':0})
mut32.loc[mut1.mutation=='I213N', 'x']=1

#import packages
import numpy
import pandas
from scipy.optimize import minimize
from scipy.stats import norm
from plotnine import *

ggplot(mut12,aes(x='x',y='y'))+geom_point()+theme_classic()

#null hypothesis likelihood equation
def fun1a(p,obs):
    B0=p[0]
    sigma=p[1]
    expected=B0
    nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
    return nll

#Alternative hypothesis likelihood equation
def fun1b(p,obs):
    B02=p[0]
    B12=p[1]
    sigma2=p[1]
    expected2=B02+B12*obs.x
    nll=-1*norm(expected2,sigma2).logpdf(obs.y).sum()
    return nll

#estimaitng parameters by minimizing the nll 
initialVals1=numpy.array([1,1,1])

fitNull=minimize(fun1a,initialVals1, method="Nelder-Mead",options={'disp': True}, args=mut12)
fitAlt=minimize(fun1b,initialVals1, method="Nelder-Mead",options={'disp': True}, args=mut12)

print(fitNull.x)
print(fitAlt.x)

from scipy.stats import chi2
D=(2*(fitAlt.fun-fitNull.fun))
mut12answer=1-chi2.cdf(x=D,df=1)
print(mut12answer)


## currently this does not produce the correct p value according to the p value sent to us from Stuart.
## the code here includes df subsets for all 3 mutations but does not find the p value for all 3. only 1
## going to correct any mistakes before adding the code for mutation 2 and mutation 3