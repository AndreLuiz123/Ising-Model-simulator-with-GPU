import matplotlib.pyplot as plt
import numpy as np

statistics_file = open("statistics.txt")
statistics_content = statistics_file.read()
statistics_array = statistics_content.split("\n")

x_axis = []
y_energy_axis = []
y_spin_axis = []
c = 0
for content in statistics_array:
    data = content.split(",")
    if len(data) > 1 and len(data)<4:
        x_axis.append(data[0])
        y_energy_axis.append(data[1])
        y_spin_axis.append(data[2])
    c = c+1

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
ax = axes[0]
ax.set_ylim([-1, 1])
ax.plot(np.asarray(x_axis, int),np.asarray(y_spin_axis, float))
ax.ticklabel_format(useOffset=False,style="plain")
ax = axes[1]
ax.plot(np.asarray(x_axis, int),np.asarray(y_energy_axis, float))
ax.ticklabel_format(useOffset=False,style="plain")
fig.tight_layout()
plt.savefig("grafico_exemplo.png", dpi=300)
plt.show()

