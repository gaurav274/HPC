P = [1,2,4,8,16] 
T1 = [0.036144, 0.0182779, 0.0100598, 0.00512886, 0.00260901]

N = [80000, 160000, 320000, 640000, 1280000, 2560000, 5120000, 10240000, 20480000]
T2 = [1.20163e-04, 2.11000e-04, 3.81947e-04, 7.00951e-04, 1.35088e-03, 2.66695e-03, 5.24998e-03, 1.04411e-02, 2.08800e-02]   

import matplotlib.pyplot as plt 
import numpy as np

fig, ax = plt.subplots() 
ax.plot(P, T1, marker = '.', linestyle = '-')
ax.set(xlabel='#Processors (P)', ylabel='Time (s)', title='P vs T [N = 5000000, C = 100]')
plt.show() 

fig, ax = plt.subplots() 
ax.plot(np.log2(N), T2, marker = '.', linestyle = '-')
ax.set(xlabel='log of #Problem Size (N)', ylabel='Time (s)', title='N vs T [P = 8, C = 100]', [])
plt.show() 