import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
from pylab import *
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
sformatter=ScalarFormatter(useOffset=True,useMathText=True)
sformatter.set_scientific(True)
sformatter.set_powerlimits((-2,3))

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(6,5))
fig = plt.figure(1)
ax = fig.add_axes([0.15,0.12,0.8,0.8])

mydata = np.loadtxt('hypercubedata.txt',skiprows=0,unpack=False)
x=mydata[0]
y=mydata[1]
mydataHO = np.loadtxt('coulomb.txt',skiprows=0,unpack=False)
xHO=mydataHO[0]
yHO=mydataHO[1]


print(x)

#x = np.arange(0,40.1,1)
#y =100.0* x**2
#z=0.5*y
#Use linestyle='None' for no line...
plt.plot(x,y,linestyle='None',color='b',markersize=8, marker='o', markerfacecolor='b', markeredgecolor='b')
plt.plot(xHO,yHO,linestyle='None',color='r',markersize=8, marker='s', markerfacecolor='r', markeredgecolor='r')

#plt.plot(x,z,linestyle=linestyles[1],linewidth=2,color='k',markersize=8, marker=markerstyles[3], markerfacecolor='r', markeredgecolor=colors[3])

#plt.semilogy(x,y)

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(-2,2,0.5), minor=False)
ax.set_xticklabels(np.arange(-2,2,0.5), minor=False, family='serif')
ax.set_xticks(np.arange(-2,2,0.1), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.xlim(-1.0,1.0)

ax.set_yticks(np.arange(-2,2,0.5), minor=False)
ax.set_yticklabels(np.arange(-2,2,0.5), minor=False, family='serif')
ax.set_yticks(np.arange(-2,2,0.1), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.ylim(-1.0,1.0)

plt.xlabel('$x$', fontsize=18, weight='normal')
plt.ylabel('$y$',fontsize=18, weight='normal')
#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('hypercube.pdf',format='pdf')
os.system('open -a Preview hypercube.pdf')
quit()
