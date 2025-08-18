import matplotlib.pyplot as plt
import numpy as np

statistics_file = open("statistics_calor_especifico.txt")
statistics_content = statistics_file.read()
statistics_array = statistics_content.split("\n")

x_axis = []
y_calor_especifico = []
c = 0
for content in statistics_array:
    data = content.split(",")
    if len(data) > 1 and len(data)<4:
        x_axis.append(data[0])
        y_calor_especifico.append(data[1])
    c = c+1

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
ax = axes[0]
ax.plot(np.asarray(x_axis,float),np.asarray(y_calor_especifico, float),'o',color='black')
ax.ticklabel_format(useOffset=False,style="plain")
fig.tight_layout()
plt.savefig("grafico_exemplo.png", dpi=300)
plt.show()

