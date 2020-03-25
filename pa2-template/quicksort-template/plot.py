P = [1,2,4,8,16] 
T1 = [8.56955, 8.94425, 5.11383, 4.02396, 2.40716] #100000000
T11 = [0.743138, 0.774818, 0.421705, 0.357514, 0.204813] #10000000

N = [10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000, 10000000000]
T222 = []
T2 = [0.000689983, 0.000705957, 0.000514984, 0.000674963, 0.00314403, 0.0286729, 0.35167, 4.14648, 40.0287, 53.0166] #P=8
T22 = [0.00104904, 0.00119209, 0.00109696, 0.00116801, 0.00250506, 0.020416, 0.204142, 2.40049, 31.6262, 28.7979] #P=16

#p=8
N2 = [1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000,9000000, 10000000, 11000000, 12000000, 13000000, 14000000, 15000000, 16000000, 17000000, 18000000, 19000000,20000000]

T3 = [0.0286288, 0.0635769, 0.0923579, 0.128511, 0.159278, 0.196721, 0.244509, 0.265889, 0.317802, 0.367397, 0.394167, 0.445023, 0.440629, 0.476861, 0.508993, 0.554705, 0.57299, 0.606249, 0.644764, 0.666871]

import matplotlib.pyplot as plt 
import numpy as np

fig, ax = plt.subplots() 
ax.plot(P, T11, marker = '.', linestyle = '-', label='P vs T [N = 10000000]')
ax.plot(P, T1, marker = '.', linestyle = '-', label='P vs T [N = 100000000]')
ax.set(xlabel='#Processors (P)', ylabel='Time (s)', title='P vs T')
leg1 = ax.legend(loc='upper right')
plt.show()

fig, ax = plt.subplots() 
ax.plot(P, T11, marker = '.', linestyle = '-', label='P vs T [N = 10000000]')
#ax.plot(P, T1, marker = '.', linestyle = '-', label='P vs T [N = 100000000]')
ax.set(xlabel='#Processors (P)', ylabel='Time (s)', title='P vs T')
leg1 = ax.legend(loc='upper right')
plt.show() 

fig, ax = plt.subplots() 
#ax.plot(P, T11, marker = '.', linestyle = '-', label='P vs T [N = 10000000]')
ax.plot(P, T1, marker = '.', linestyle = '-', label='P vs T [N = 100000000]')
ax.set(xlabel='#Processors (P)', ylabel='Time (s)', title='P vs T')
leg1 = ax.legend(loc='upper right')
plt.show()

#############################################

fig, ax = plt.subplots() 
ax.plot(N, T2, marker='.', linestyle='-', label='N vs T [P = 8]')
#ax.plot(N, np.log(T22), marker = '.', linestyle = '-', label='N vs T [P = 16]')
ax.set(xlabel='#Problem Size (N)', ylabel='Time (s)', title='N vs T')
leg1 = ax.legend(loc='upper right')
plt.show()

fig, ax = plt.subplots() 
ax.plot(N2, T3, marker='.', linestyle='-', label='N vs T [P = 8]')
#ax.plot(N, np.log(T22), marker = '.', linestyle = '-', label='N vs T [P = 16]')
ax.set(xlabel='#Problem Size (N)', ylabel='Time (s)', title='N vs T')
leg1 = ax.legend(loc='upper right')
plt.show() 

#fig, ax = plt.subplots() 
#ax.plot(np.log2(N), T2, marker = '.', linestyle = '-')
#ax.set(xlabel='log of #Problem Size (N)', ylabel='Time (s)', label='log(N) vs T [P = 8, C = 100]')
#plt.show() 

# fig, ax = plt.subplots() 
# ax.plot(np.log2(N), T2, marker = '.', linestyle = '-', label='log(N) vs T [P = 8, C = 100]')
# ax.plot(np.log2(N), T3, marker = '.', linestyle = '-', label='log(N) vs T [P = 16, C = 100]')
# ax.set(xlabel='log of #Problem Size (N)', ylabel='Time (s)', title='log(N) vs T')
# leg1 = ax.legend(loc='upper left')
# plt.show() 
