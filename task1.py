from qutip import *
import matplotlib.pyplot as plt
import numpy as np

import time
#import psutil

class system:
    def __init__(self, ham_matrix, psi_init, fock_dim):
        self.a = []
        for i in range(0, len(ham_matrix)):
            op_list = []
            for j in range(0, len(ham_matrix)):
                if j == i:
                    op_list.append(destroy(fock_dim))
                else:
                    op_list.append(qeye(fock_dim))
            self.a.append(tensor(op_list))

        for i in range(0, len(ham_matrix)):
            for j in range(0, len(ham_matrix)):
                if i == 0 and j == 0:
                    self.H = ham_matrix[0][0] * self.a[0].dag() * self.a[0]
                else:
                    self.H = self.H + ham_matrix[i][j] * self.a[i].dag() * self.a[j]

    def gen_c_ops(self, gamma):
        c_ops_list = []
        for i in range(0, len(self.a)):
            c_ops_list.append(gamma * self.a[i])
        return c_ops_list

    def gen_cormatrix_2(self, psi):
        cormatrix = np.zeros((len(self.a), len(self.a)), dtype=complex)
        for i in range(0, len(self.a)):
            for j in range(0, len(self.a)):
                temp_op = self.a[i].dag() * self.a[j].dag() * self.a[i] * self.a[j]
                cormatrix[i][j] = temp_op.matrix_element(psi.dag(), psi)
        return cormatrix

    def gen_cormatrix_2_loss(self, psi):
        cormatrix = np.zeros((len(self.a), len(self.a)), dtype=complex)
        for i in range(0, len(self.a)):
            for j in range(0, len(self.a)):
                temp_op = self.a[i].dag() * self.a[j].dag() * self.a[i] * self.a[j]
                temp_op_j = psi * temp_op
                cormatrix[i][j] = temp_op_j.tr()
        return cormatrix

def psi_init(psi_vector, fock_dim):
    psi_list = []
    for i in range(0,len(psi_vector)):
        psi_list.append(fock(fock_dim,psi_vector[i]))
    print("LOOOOOOOOOL")
    print(psi_list)
    return tensor(psi_list)

def uniform_array(num_wg, b, k):
    ham_matrix = np.zeros((num_wg,num_wg))
    for i in range(0,num_wg):
        for j in range(0,num_wg):
            if i == j:
                ham_matrix[i][j] = b
            elif i == j+1 or j == i+1:
                ham_matrix[i][j] = k
    return ham_matrix

n = 6
psi_vector = np.zeros(n)
psi_vector = psi_vector.astype(int)
psi_vector[2] = 1
psi_vector[3] = 1

psi0 = psi_init(psi_vector, 3)
print("________PSI")
#print(psi0)
print("________PSI0")

hm = uniform_array(n, 1, 1)
print(hm)
print("________HM")
qws = system(hm, psi0, 3)
print("________ASTART")
#print(qws.a)
print("________A")
times = np.linspace(0,5,100)

start = time.time()
#result_un = sesolve(qws.H, psi0, times, [])
end = time.time()
print(end - start)

#result_loss = mesolve(qws.H, psi0, times, qws.gen_c_ops(0.25), [])


'''n = list(range(2,11))
calc_time = []
mem_use = []
mem_bg = 0
process = psutil.Process(os.getpid())
for i in range(0,len(n)):
    psi_vector = np.zeros(n[i])
    psi_vector = psi_vector.astype(int)
    psi_vector[int(round(n[i]/2))] = 1
    psi_vector[int(round(n[i]/2-1))] = 1
    psi0 = psi_init(psi_vector, 3)
    hm = uniform_array(n[i], 1, 1)
    qws = system(hm, psi0, 3)
    times = np.linspace(0,5,100)
    start = time.time()
    result_un = sesolve(qws.H, psi0, times, [])
    end = time.time()
    calc_time.append(end - start)
    mem_use.append(process.memory_info().rss/1024/1024)
    print(process.memory_info().rss/1024/1024)
    print(end - start)
    
plt.semilogy(n, calc_time, label = 'Two-photons')
plt.legend()
plt.ylabel('Calculation time (sec)')
plt.xlabel('Number of modes')
plt.show()
'''