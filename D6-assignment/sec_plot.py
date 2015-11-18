import numpy as np
import matplotlib.pyplot as plt

mat = np.loadtxt("transp_res.dat")

th = mat[:,0]
t = mat[:,1]

# th = np.array(nth, dtype = int)
# sp = np.array(t, dtype = float)
speed = np.array([t[0] / t[i] for i in range(len(t))])

print(th)
print(speed)

plt.figure()
plt.plot(th, speed, '-o', label = 'mes')
plt.plot(th, th, label = 'att')
plt.xlabel("N_Threads")
plt.ylabel("time (s)")
plt.title("Transpose scaling")
plt.xlim(th.min(), th.max())

plt.legend(bbox_to_anchor=(.22,1))
# plt.show()
plt.savefig("transp_scaling.png")
