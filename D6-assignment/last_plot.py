import numpy as np
import matplotlib.pyplot as plt

mat  = np.loadtxt("data_scaling_omp_th_zero")
mat1 = np.loadtxt("data_scaling_omp_th_two")
mat2 = np.loadtxt("data_scaling_omp_th_five")

th = mat[:,0]
t1 = mat[:,1]
t2 = mat1[:,1]
t3 = mat2[:,1]

# th = np.array(nth, dtype = int)
# sp = np.array(t, dtype = float)
speed  = np.array([t1[0] / t1[i] for i in range(len(t1))])
speed1 = np.array([t2[0] / t2[i] for i in range(len(t2))])
speed2 = np.array([t3[0] / t3[i] for i in range(len(t3))])

print(th)
print(speed)

plt.figure()
plt.plot(th, speed,  '-o', label = 'OMP_NUM_THREADS = 0')
plt.plot(th, speed1, '-o', label = 'OMP_NUM_THREADS = 2')
plt.plot(th, speed2, '-o', label = 'OMP_NUM_THREADS = 5')
plt.plot(th, th, label = 'att')
plt.xlabel("N_Threads")
plt.ylabel("time (s)")
plt.title("Transpose scaling")
plt.xlim(th.min(), th.max())

plt.legend(bbox_to_anchor=(.55, 1))
# plt.show()
plt.savefig("mpi_vs_omp_scaling.png")
