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


x_axis_a0, y_calor_especifico_a0 = plot_shart("a_0.0_M_176_N_10.txt",1.3,3)
x_axis_a15, y_calor_especifico_a15 = plot_shart("a_0.15_M_176_N_10.txt",1.0,3)
x_axis_a2, y_calor_especifico_a2 = plot_shart("a_0.20_M_176_N_10.txt",0.9,2.7)
x_axis_a25, y_calor_especifico_a25 = plot_shart("a_0.25_M_176_N_10.txt",0.5,3)
x_axis_a285, y_calor_especifico_a285 = plot_shart("a_0.285_M_176_N_10.txt",0,3)

x_axis_a376, y_calor_especifico_a376 = plot_shart("a_0.376_M_176_N_10.txt",0.3,2.5)
x_axis_a637, y_calor_especifico_a637 = plot_shart("a_0.637_M_176_N_10.txt",0.9,2.9)
x_axis_a8, y_calor_especifico_a8 = plot_shart("a_0.800_M_352_N_10.txt",1.3,3)
x_axis_a1, y_calor_especifico_a1 = plot_shart("a_1_M_176_N_10.txt",1.3,3)




fig, axes = plt.subplots(1,2, figsize=(8, 4))
ax = axes[0]
ax.plot(np.asarray(x_axis_a0,float),np.asarray(y_calor_especifico_a0, float),'-',label="a = 0.0",color='black')
ax.plot(np.asarray(x_axis_a15,float),np.asarray(y_calor_especifico_a15, float),'-',label="a = 0.15",color='blue')
ax.plot(np.asarray(x_axis_a2,float),np.asarray(y_calor_especifico_a2, float),'-',label="a = 0.2",color='red')
ax.plot(np.asarray(x_axis_a25,float),np.asarray(y_calor_especifico_a25, float),'-',label="a = 0.25",color='orange')
ax.plot(np.asarray(x_axis_a285,float),np.asarray(y_calor_especifico_a285, float),'-',label="a = 0.285",color='green')
ax.ticklabel_format(useOffset=False,style="plain")
ax.set_xlabel("Temperatura")
ax.set_ylabel("Calor específico")
ax.legend()
ax = axes[1]
ax.plot(np.asarray(x_axis_a376,float),np.asarray(y_calor_especifico_a376, float),'-',label = "a = 0.376",color='purple')
ax.plot(np.asarray(x_axis_a637,float),np.asarray(y_calor_especifico_a637, float),'-',label = "a = 0.637",color='brown')
ax.plot(np.asarray(x_axis_a8,float),np.asarray(y_calor_especifico_a8, float),'-',label = "a = 0.852",color='gold')
ax.plot(np.asarray(x_axis_a1,float),np.asarray(y_calor_especifico_a1, float),'-',label = "a = 1.0",color='grey')
ax.ticklabel_format(useOffset=False,style="plain")
ax.set_xlabel("Temperatura")
ax.set_ylabel("Calor específico")
ax.legend()
fig.tight_layout()
plt.savefig("grafico_exemplo.png", dpi=300)
plt.show()

