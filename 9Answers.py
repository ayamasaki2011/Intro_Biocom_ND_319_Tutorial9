#Exercise 09
###################################################################

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

#plot the values for WT at 0 and mutation 1,2,3 at 1 --- dont have to plot anything but it helped me visualize
ggplot(mut12,aes(x='x',y='y'))+geom_point()+theme_classic()
ggplot(mut22,aes(x='x',y='y'))+geom_point()+theme_classic()
ggplot(mut32,aes(x='x',y='y'))+geom_point()+theme_classic()

#null hypothesis likelihood equation
def null(p,obs):
    B0=p[0]
    sigma=p[1]
    expected=B0
    nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
    return nll

#Alternative hypothesis likelihood equation
def alt(p,obs):
    B0=p[0]
    B1=p[1]
    sigma=p[2]
    expected=B0+B1*obs.x
    nll=-1*norm(expected,sigma).logpdf(obs.y).sum()
    return nll

#estimaitng parameters by minimizing the nll 
initialVals1=numpy.array([1,1,1])

fitNull=minimize(null,initialVals1, method="Nelder-Mead",options={'disp': True}, args=mut12)
fitAlt=minimize(alt,initialVals1, method="Nelder-Mead",options={'disp': True}, args=mut12)

print(fitNull.x)
print(fitAlt.x)

from scipy.stats import chi2
D=(2*(fitNull.fun-fitAlt.fun))
mut12answer=(1-chi2.cdf(x=D,df=1))
print('mutation M124K p value')
print(mut12answer)


## for mutation 2
fitNull=minimize(null,initialVals1, method="Nelder-Mead",options={'disp': True}, args=mut22)
fitAlt=minimize(alt,initialVals1, method="Nelder-Mead",options={'disp': True}, args=mut22)

print(fitNull.x)
print(fitAlt.x)

from scipy.stats import chi2
D=(2*(fitNull.fun-fitAlt.fun))
mut22answer=1-chi2.cdf(x=D,df=1)
print('mutation V456D p value')
print(mut22answer)

##for mutation 3

fitNull=minimize(null,initialVals1, method="Nelder-Mead",options={'disp': True}, args=mut32)
fitAlt=minimize(alt,initialVals1, method="Nelder-Mead",options={'disp': True}, args=mut32)

print(fitNull.x)
print(fitAlt.x)

from scipy.stats import chi2
D=(2*(fitNull.fun-fitAlt.fun))
mut32answer=1-chi2.cdf(x=D,df=1)
print('mutation I213N p value')
print(mut32answer)


############################################################################################################
# question number 2

#import and read file for question 2
import pandas
file2=pandas.read_csv("MmarinumGrowth.csv",header=0,sep=",")

#again, dont need to plot but it helped me with the dataset
ggplot(file2,aes(x='S',y='u'))+geom_point()+theme_classic()

#changed the typical equation to the equation given -see expected=...
## and aligned the correct variable names with the given symbols from the exercise question

def nllike(p,obs):
    umax=p[0]
    Ks=p[1]
    sigma=p[2]
    expected=umax*obs.S/(obs.S+Ks)
    nll=-1*norm(expected,sigma).logpdf(obs.u).sum()
    return nll

#estimaitng parameters by minimizing the nll 
guess2=numpy.array([1, 1, 1])

fitNull2=minimize(nllike,guess2, method="Nelder-Mead",options={'disp': True}, args=file2)

print(fitNull2.x)
print('umax Ks  sigma')

###################################################################################################################

#question 3

#import and read the file for question 3
import pandas
file3=pandas.read_csv("leafDecomp.csv",header=0,sep=",")

#same notes as before - dont have to do a plot but it helped me
ggplot(file3,aes(x='Ms',y='decomp'))+geom_point()+theme_classic()

#constant fit equation
def constantfit(p,obs):
    a=p[0]
    sigma=p[2]
    expected=a
    nll=-1*norm(expected,sigma).logpdf(obs.decomp).sum()
    return nll

guess3=numpy.array([1, 1, 1])

fitNull3constant=minimize(constantfit,guess3, method="Nelder-Mead",options={'disp': True}, args=file3)

print('constant fit')
print(fitNull3constant.x)



### linear fit
def linearfit(p,obs):
    a=p[0]
    b=p[1]
    sigma=p[2]
    expected=a+b*obs.Ms 
    nll=-1*norm(expected,sigma).logpdf(obs.decomp).sum()
    return nll

guess3=numpy.array([1, 1, 1])

fitNull3linear=minimize(linearfit,guess3, method="Nelder-Mead",options={'disp': True}, args=file3)

print('linear fit')
print(fitNull3linear)

from scipy.stats import chi2
D=(2*(fitNull3constant.fun-fitNull3linear.fun))
linearVSconstant=(1-chi2.cdf(x=D,df=1))
print('linear vs constant p value')
print(linearVSconstant)

### quad fit
def quadfit(p,obs):
    a=p[0]
    b=p[1]
    c=p[2]
    sigma=p[3]
    expected=a+b*obs.Ms+c*((obs.Ms)*(obs.Ms))
    nll=-1*norm(expected,sigma).logpdf(obs.decomp).sum()
    return nll

guess3=numpy.array([180, 15.7, -0.11, 10])

fitNull3quad=minimize(quadfit,guess3, method="Nelder-Mead",options={'disp': True}, args=file3)

print('quadratic fit')
print(fitNull3quad)

from scipy.stats import chi2
D=(2*(fitNull3constant.fun-fitNull3quad.fun))
quadVSconstant=(1-chi2.cdf(x=D,df=1))
print('quadratic vs constant p value')
print(quadVSconstant)
