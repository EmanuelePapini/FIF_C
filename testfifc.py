"""
python file to test the output of TestFif.c

first compile the code with "make"

then execute ./test

then run python3 ./testfifc.py

"""

import numpy as np
import pylab as plt

numIMF = 3
plt.figure()
sig = np.loadtxt('signal.dat')

for i in range(1,numIMF+1):
    imf = np.loadtxt('IMF'+str(i)+'.dat')
    plt.plot(imf,label='IMF '+str(i))

res = np.loadtxt('residual.dat')
plt.plot(res,label='residual')

plt.plot(sig,label='original signal')
plt.legend()
plt.show()
