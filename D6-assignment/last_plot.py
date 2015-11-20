import numpy as np
import matplotlib.pyplot as plt

mat  = np.loadtxt("data_scaling_omp_th_zero")
mat1 = np.loadtxt("data_scaling_omp_th_two")
mat2 = np.loadtxt("data_scaling_omp_th_five")
mat3 = np.loadtxt("data_scaling_omp_th_ten")

th = mat[:,0]
t1 = mat[:,1]
t2 = mat1[:,1]
t3 = mat2[:,1]
t4 = mat3[:,1]

speed  = np.array([t1[0] / t1[i] for i in range(len(t1))])
speed1 = np.array([t1[0] / t2[i] for i in range(len(t2))])
speed2 = np.array([t1[0] / t3[i] for i in range(len(t3))])
speed3 = np.array([t1[0] / t4[i] for i in range(len(t4))])

print(th)
print(t1)
print(t2)

n_mes = 6
bar_width = .2
index = np.arange(1, n_mes + 1);
print(index)

plt.figure()

plt.bar(index, t1, bar_width, label = 'OMP_NUM_THREADS = 0')
plt.bar(index + bar_width, t2, bar_width, color = 'r', label = 'OMP_NUM_THREADS = 2')
plt.bar(index + 2 * bar_width, t3, bar_width, color = 'g', label = 'OMP_NUM_THREADS = 5')
plt.bar(index + 3 * bar_width, t4, bar_width, color = 'c', label = 'OMP_NUM_THREADS = 10')

plt.title("MPI vs OpenMP timing")
plt.xlabel("Number of Nodes")
plt.ylabel("Time (s)")
plt.xticks(index + 2 * bar_width, index)
plt.legend()
plt.savefig("mpi_vs_omp_timing.png")
plt.close('all')

plt.figure()
plt.plot(th, speed,  '-o', label = 'OMP_NUM_THREADS = 0')
plt.plot(th, speed1, '-o', label = 'OMP_NUM_THREADS = 2')
plt.plot(th, speed2, '-o', label = 'OMP_NUM_THREADS = 5')
plt.plot(th, speed3, '-o', label = 'OMP_NUM_THREADS = 10')

plt.xlabel("Number of nodes")
plt.ylabel("Speedup")
plt.title("MPI vs OpenMP scaling")
plt.xlim(th.min(), th.max())

plt.legend(bbox_to_anchor=(.575, 1))
# plt.show()
plt.savefig("mpi_vs_omp_scaling.png")
