import matplotlib.pyplot as plt
import numpy as np


def occupation_derivative(a,temp):
    """
                   df
        Computes ------ where f = Bose occupation factor
                   dT
    """
    beta = 100/(8.617333262*temp)
    return np.exp(a*beta)/((np.exp(a*beta)-1)**2) *a*beta/temp


bose = lambda e,t : 1/(np.exp(100*e/(8.617333262*t))-1)


en = np.linspace(1,250,1000)
Temp = [10, 15, 20, 50, 100, 200, 300] #K

data = np.zeros((len(en),len(Temp),2),dtype=np.float64)

for i in range(len(en)):
    for j in range(len(Temp)):
        data[i,j,0] = occupation_derivative(en[i],Temp[j])
        data[i,j,1] = bose(en[i],Temp[j])

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

for i,t in enumerate(Temp):
    ax1.plot(en, data[:,i,0], linestyle='-.' , label=r"$\frac{\partial f}{\partial T}_{%3d K}$"%(t))
    ax2.plot(en, data[:,i,1][::-1],label=r"$%3d K$ Bose"%(t))
ax1.set_ylabel("Derivatives")
ax2.set_ylabel("Bose")
ax1.legend(loc=8)
ax2.legend(loc=9)
plt.show()
