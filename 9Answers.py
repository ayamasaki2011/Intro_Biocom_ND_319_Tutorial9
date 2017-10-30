#Exercise 09
###################################################################3

# Question 1 zebrafish mutations

#import pandas and read csv file
import pandas
file=pandas.read_csv("ponzr1.csv",header=0,sep=",")

#subset each mutation values with WT values
mut1=file.loc[file.mutation.isin(['WT', 'M124K']),:]
mut2=file.loc[file.mutation.isin(['WT', 'V456D']),:]
mut3=file.loc[file.mutation.isin(['WT', 'I213N']),:]

#mutation 1 subset into a new dataframe and changes the x column to 0's and 1's
mut12=pandas.DataFrame({'y':mut1.ponzr1Counts, 'x':0})
mut12.loc[mut1.mutation=='M124K', 'x']=1

#mutation 2 subset into a new dataframe and changes the x column to 0's and 1's
mut22=pandas.DataFrame({'y':mut2.ponzr1Counts, 'x':0})
mut22.loc[mut2.mutation=='V456D', 'x']=1

#mutation 3 subset into a new dataframe and changes the x column to 0's and 1's
mut32=pandas.DataFrame({'y':mut3.ponzr1Counts, 'x':0})
mut32.loc[mut3.mutation=='I213N', 'x']=1

#import packages
import numpy
import pandas
from scipy.optimize import minimize
from scipy.stats import norm
from plotnine import *

#plot the values for WT at 0 and mutation 1 at 1 --- not sure if we need to do this
ggplot(mut12,aes(x='x',y='y'))+geom_point()+theme_classic()

#null hypothesis likelihood equation
def nllike(p,obs):
    B0=p[0]
    sigma=p[1]
    expected=B0
    nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
    return nll

#Alternative hypothesis likelihood equation
def nllike2(p,obs):
    B0=p[0]
    B1=p[1]
    sigma=p[2]
    expected=B0+B1*obs.x
    nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
    return nll

#estimaitng parameters by minimizing the nll 
initialVals1=numpy.array([1,1,1])

fitNull=minimize(nllike,initialVals1, method="Nelder-Mead",options={'disp': True}, args=mut12)
fitAlt=minimize(nllike2,initialVals1, method="Nelder-Mead",options={'disp': True}, args=mut12)

print(fitNull.x)
print(fitAlt.x)

from scipy.stats import chi2
D=(2*(fitAlt.fun-fitNull.fun))
mut12answer=(1-chi2.cdf(x=D,df=1))
print('mutation M124K p value')
print(mut12answer)


## for mutation 2
fitNull=minimize(fun1a,initialVals1, method="Nelder-Mead",options={'disp': True}, args=mut22)
fitAlt=minimize(fun1b,initialVals1, method="Nelder-Mead",options={'disp': True}, args=mut22)

print(fitNull.x)
print(fitAlt.x)

from scipy.stats import chi2
D=(2*(fitAlt.fun-fitNull.fun))
mut22answer=1-chi2.cdf(x=D,df=1)
print('mutation V456D p value')
print(mut22answer)

##for mutation 3

fitNull=minimize(fun1a,initialVals1, method="Nelder-Mead",options={'disp': True}, args=mut32)
fitAlt=minimize(fun1b,initialVals1, method="Nelder-Mead",options={'disp': True}, args=mut32)

print(fitNull.x)
print(fitAlt.x)

from scipy.stats import chi2
D=(2*(fitAlt.fun-fitNull.fun))
mut32answer=1-chi2.cdf(x=D,df=1)
print('mutation I213N p value')
print(mut32answer)


## currently this does not produce the correct p value according to the p value sent to us from Stuart.
## the code here includes df subsets for all 3 mutations and does find the p value for all 3. 
## its just not right

############################################################################################################
# question number 2
import pandas
file2=pandas.read_csv("MmarinumGrowth.csv",header=0,sep=",")

ggplot(file2,aes(x='S',y='u'))+geom_point()+theme_classic()

def nllike(p,obs):
    B0=p[0]
    B1=p[1]
    sigma=p[2]
    expected=B0
    nll=-1*norm(expected,sigma).logpdf(obs.S).sum()
    return nll

guess2=numpy.array([1000, 1000, 1000])

fitNull2=minimize(nllike,guess2, method="Nelder-Mead",options={'disp': True}, args=file2)

print(fitNull2.x)

###################################################################################################################

#question 3
import pandas
file3=pandas.read_csv("leafDecomp.csv",header=0,sep=",")

ggplot(file3,aes(x='Ms',y='decomp'))+geom_point()+theme_classic()
#constant fit 
def nllike(p,obs):
    B0=p[0]
    B1=p[1]
    sigma=p[2]
    expected=B0
    nll=-1*norm(expected,sigma).logpdf(obs.decomp).sum()
    return nll

guess3=numpy.array([1, 1, 1])

fitNull3constant=minimize(nllike,guess3, method="Nelder-Mead",options={'disp': True}, args=file3)

print('constant fit')
print(fitNull3constant.x)

### linear fit
def nllike(p,obs):
    B0=p[0]
    B1=p[1]
    sigma=p[2]
    expected=B0+B1*obs.Ms 
    nll=-1*norm(expected,sigma).logpdf(obs.decomp).sum()
    return nll

fitNull3linear=minimize(nllike,guess3, method="Nelder-Mead",options={'disp': True}, args=file3)

print('linear fit')
print(fitNull3linear)

### quad fit --- doesn't work
def nllike(p,obs):
    B0=p[0]
    B1=p[1]
    B2=p[2]
    sigma=p[2]
    expected=B0+B1*obs.Ms+B2*((obs.Ms)*(obs.Ms))
    nll=-1*norm(expected,sigma).logpdf(obs.decomp).sum()
    return nll

fitNull3quad=minimize(nllike,guess3, method="Nelder-Mead",options={'disp': True}, args=file3)

print('quadratic fit')
print(fitNull3quad)

