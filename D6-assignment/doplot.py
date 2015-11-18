import numpy as np
import matplotlib.pyplot as plt

# fl = open("toparse")
fl = open("new_parse")

nth = []
t = []

for line in fl:
    if len(line) > 1 and line[1] == "T":
        a = line.split()
        # print(line)
        # print(a[2])
        t.append(a[2])
    elif len(line) > 1 and line[0] == "\t" and line[1] == "N":
        # print(line)
        a = line.split()
        # print(a[-1])
        nth.append(a[-1])

th = np.array(nth, dtype = int)
sp = np.array(t, dtype = float)
speed = np.array([sp[0] / sp[i] for i in range(len(sp))])

print(th)
print(speed)

plt.figure()
plt.plot(th, speed, '-o', label = 'mes')
plt.plot(th, th, label = 'att')
plt.xlabel("N_Threads")
plt.ylabel("time (s)")
plt.title("Heat equation scaling")
plt.xlim(th.min(), th.max())

# print(th.min(), "\t", th.max())

plt.legend(bbox_to_anchor=(.22,1))
# plt.show()
plt.savefig("heat_scaling.png")
