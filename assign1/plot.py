P = [1,2,4,8,16] 
T1 = [6.18532e-02, 3.17831e-02, 1.70200e-02, 9.57894e-03, 4.59099e-03]
N = [1280000, 2560000, 5120000, 10240000, 20480000, 40960000]

T2 = [2.53701e-03, 5.39708e-03, 9.80091e-03, 1.82111e-02, 3.56920e-02, 7.06e-02]


T3 = [
1.31083e-03,
2.57799e-03,
4.72811e-03,
9.48360e-03,
1.89359e-02,
3.68e-02]

import matplotlib.pyplot as plt 
import numpy as np

fig, ax = plt.subplots() 
ax.plot(P, T1, marker = '.', linestyle = '-')
ax.set(xlabel='#Processors (P)', ylabel='Time (s)', title='P vs T [N = 5000000, C = 100]')
plt.show() 

fig, ax = plt.subplots() 
ax.plot(N, T2, marker = '.', linestyle = '-')
ax.set(xlabel='#Problem Size (N)', ylabel='Time (s)', title='N vs T [P = 8, C = 100]')
plt.show() 

#fig, ax = plt.subplots() 
#ax.plot(np.log2(N), T2, marker = '.', linestyle = '-')
#ax.set(xlabel='log of #Problem Size (N)', ylabel='Time (s)', label='log(N) vs T [P = 8, C = 100]')
#plt.show() 

fig, ax = plt.subplots() 
ax.plot(np.log2(N), T2, marker = '.', linestyle = '-', label='log(N) vs T [P = 8, C = 100]')
ax.plot(np.log2(N), T3, marker = '.', linestyle = '-', label='log(N) vs T [P = 16, C = 100]')
ax.set(xlabel='log of #Problem Size (N)', ylabel='Time (s)', title='log(N) vs T')
leg1 = ax.legend(loc='upper left')
plt.show() 
