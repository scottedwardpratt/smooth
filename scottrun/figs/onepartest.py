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

Lambda=float(input("Enter Lambda: "))
TrainRank=int(input("Enter NThetaPoints: "))-1

print('Lambda=',Lambda,' TrainRank=',TrainRank)

for ireal in range(0,10):
  filename='../testresults/NPars1_Lambda'+str(Lambda)+'_TrainRank'+str(TrainRank)+'_real'+str(ireal)+'.txt'
  mydata = np.loadtxt(filename,skiprows=0,unpack=True)
  print('filename=',filename)
  x=mydata[0]
  yreal=mydata[1]
  y1=mydata[2]
  y2=mydata[3]
  y3=mydata[4]
  y4=mydata[5]
  y5=mydata[6]
  y6=mydata[7]
  y7=mydata[8]
  y8=mydata[9]
  y9=mydata[10]
  y10=mydata[11]

  if TrainRank==1:
    xpoint=[-0.675,0.675]
    ypoint=[yreal[6],yreal[33]]
  if TrainRank==2:
    xpoint=[-0.825,0,0.825]
    ypoint=[yreal[3],0.5*yreal[19]+0.5*yreal[20],yreal[36]]
  else:
    print('TrainRank must be 1 or 2')

  plt.plot(xpoint,ypoint,linestyle='None',color='k',markersize=8, marker='o', markerfacecolor='k', markeredgecolor='k')
  plt.plot(x,y1,linestyle='-',linewidth=1,color='r')
  plt.plot(x,y2,linestyle='-',linewidth=1,color='r')
  plt.plot(x,y3,linestyle='-',linewidth=1,color='r')
  plt.plot(x,y4,linestyle='-',linewidth=1,color='r')
  plt.plot(x,y5,linestyle='-',linewidth=1,color='r')
  plt.plot(x,y6,linestyle='-',linewidth=1,color='r')
  plt.plot(x,y7,linestyle='-',linewidth=1,color='r')
  plt.plot(x,y8,linestyle='-',linewidth=1,color='r')
  plt.plot(x,y9,linestyle='-',linewidth=1,color='r')
  plt.plot(x,y10,linestyle='-',linewidth=1,color='r')
  plt.plot(x,yreal,linestyle='-',linewidth=1,color='k')

#plt.plot(x,z,linestyle=linestyles[1],linewidth=2,color='k',markersize=8, marker=markerstyles[3], markerfacecolor='r', markeredgecolor=colors[3])

#plt.semilogy(x,y)

ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(-2,2,0.5), minor=False)
ax.set_xticklabels(np.arange(-2,2,0.5), minor=False, family='serif')
ax.set_xticks(np.arange(-2,2,0.1), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.xlim(-1.0,1.0)

title("$\Lambda=$"+str(Lambda))

#ax.set_yticks(np.arange(-10,10,1.0), minor=False)
#ax.set_yticklabels(np.arange(-10,10,1.0), minor=False, family='serif')
#ax.set_yticks(np.arange(-10,10,0.5), minor=True)
#plt.ylim(-3,3)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$\\theta$', fontsize=18, weight='normal')
plt.ylabel('Y',fontsize=18)
#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('onepartest_Lambda'+str(Lambda)+'.pdf',format='pdf')
os.system('open -a Preview onepartest_Lambda'+str(Lambda)+'.pdf')
#plt.show()
quit()
