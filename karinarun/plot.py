import matplotlib
import matplotlib.pyplot as plt
import numpy as np

filename = "emulator.txt"
mydata = np.genfromtxt(filename, delimiter = None)[:,:]
#plt.figure()
x=np.linspace(0,10,1000)
y = np.squeeze(np.cos(0.1*x))
#plt.subplot(211)
plt.plot(x,y, label=r"$f(x) = \cos(0.1x)$", linestyle="dotted")
plt.plot(x,mydata, 'g^', label='emulator')
#plt.plot(mydata1[:22,0],mydata1[:22,3], 'bs', label='CFlong')
plt.legend()
plt.xlabel("$x$")
plt.ylabel("$f(x)$")
plt.title('real function')
plt.show()
