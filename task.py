import numpy as np
import matplotlib.pyplot as plt
import datetime

from qutip import *


Sx = sigmax()
Sy = sigmay()
Sz = sigmaz()
Sm = sigmam()

a = destroy(2)
b = create(2)

#H = (0.+1.j)*np.pi*(a*b)
#H1 = 2 * np.pi * 0.1 * sigmax()

phi = tensor(basis(2,0),basis(2,1))
H = tensor(b,a.dag())

print(H)
print(phi)
print(H*phi)
print(b)
print(a)
'''
print(a*phi)
print(b*phi)
print(b*a*phi)'''

#plt.show()

'''print(tensor(basis(2,0), basis(2,1)))
print(tensor(basis(2,1), basis(2,0)))
print(tensor(basis(2,1), basis(2,1)))'''
