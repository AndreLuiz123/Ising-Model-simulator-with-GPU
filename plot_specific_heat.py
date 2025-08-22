import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) == 1:
    print("Este script só funciona com um arquivo de entrada. Por exemplo: python plot_monte_carlo.py caminho/para/o/seu_arquivo.txt") 
else:
    statistics_file = open(sys.argv[1])
    statistics_content = statistics_file.read()
    statistics_array = statistics_content.split("\n")

    x_temperatura = []
    y_calor_especifico = []
    c = 0
    for content in statistics_array:
        data = content.split(",")
        if len(data) > 1 and len(data)<4:
            x_axis.append(data[0])
            y_calor_especifico.append(data[1])
        c = c+1

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(np.asarray(x_axis,float),np.asarray(y_calor_especifico, float),'o',color='black')
    ax.ticklabel_format(useOffset=False,style="plain")
    fig.tight_layout()
    plt.savefig("grafico_exemplo.png", dpi=300)
    plt.show()