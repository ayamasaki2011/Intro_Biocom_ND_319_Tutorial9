import pandas
import scipy
import scipy.integrate as spint
from plotnine import *

def CancerSim(y,t,rN,rT,Kn,Kt,aNT,aTN):
    Nn=y[0]
    Nt=y[1]
    
    dNndt=rN*(1-((Nn+aNT*Nt)/Kn))*Nn
    dNtdt=rT*(1-((Nt+aTN*Nn)/Kt))*Nt
    
    return [dNndt, dNtdt]
    
params2=(0.5,0.5,10,10,0.5,2)
params3=(0.5,0.5,10,10,0.5,0.5)
params4=(0.5,0.5,10,10,2,0.5)

N0=[0.01,0.01]
times=range(0,500)

modelSim2=spint.odeint(func=CancerSim,y0=N0,t=times,args=params2)

modelOutput2=pandas.DataFrame({"t":times,"Nn":modelSim2[:,0],"Nt":modelSim2[:,1]})

g2=ggplot(modelOutput2, aes(x="t"))
g2=g+geom_line(aes(y="Nn"))
g2=g+geom_line(aes(y="Nt"))

g2