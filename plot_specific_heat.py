import matplotlib.pyplot as plt
import numpy as np

def plot_shart(file, min_t, max_t):
    statistics_file = open(file)
    statistics_content = statistics_file.read()
    statistics_array = statistics_content.split("\n")

    x_axis = []
    y_calor_especifico = []
    c = 0
    for content in statistics_array:
        data = content.split(",")
        if len(data) > 1 and len(data)<4:
            if float(data[0]) >= min_t and float(data[0]) <= max_t:
                x_axis.append(data[0])
                y_calor_especifico.append(data[1])
        c = c+1
    statistics_file.close()
    return x_axis, y_calor_especifico


x_axis,y_calor_especifico = plot_shart("specific_heat_M_176_N_10_B_1_mcs_50000____03_09_2025_13_35_05.txt",0,3)



fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(np.asarray(x_axis,float),np.asarray(y_calor_especifico, float),'o',color='red', markersize=1)
ax.set_ylim(0.0,1.5)
#intervalo_x = np.arange(0.6, 3.1, 0.3)
ax.set_xlim(0.6,3.0)
tick_locations = np.arange(0.6, 3.1, 0.3)
ax.set_xticks(tick_locations)
tick_locations = np.arange(0.0, 1.5, 0.3)
ax.set_yticks(tick_locations)

ax.ticklabel_format(useOffset=False,style="plain")
fig.tight_layout()
plt.savefig("grafico_exemplo2.png", dpi=300)
plt.show()

