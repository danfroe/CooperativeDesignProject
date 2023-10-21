# import module
import numpy as np
import matplotlib.pyplot as plt


# initialize x and y coordinates
x = [3,5,10,15,20,25,30,35,40] #adiabatic operation
y = [42,19.5,6.9,3.8,2.5,1.7,1.3,1.1,0.9] #adiabatic operation

xiso = [3,5,10,15,20,25,30,35,40] #isothermal operation
yiso = [171.8,79.8,28.2,15.3,10,7.1,5.5,4.3,3.6] #isothermal operation


# plot scatter plot with x and y data
plt.scatter(x, y, label='Adiabatic', facecolors='none', edgecolors='royalblue', marker='s',  s = 140)
plt.scatter(25, 1.7,facecolors='royalblue', edgecolors='royalblue', marker='s',  s = 140)
plt.scatter(x, yiso,label='Isothermal', facecolors='none', edgecolors='crimson', marker='o',  s = 140)

plt.yticks([0,5,10,15,20,25])
plt.xticks([0,10,20,30,40,50])

plt.ylim(0,25)
plt.xlim(0,50)

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.rcParams.update({'font.size': 28})
plt.ylabel('$V_{\mathrm{R}}$ / ($\mathrm{m^3})$')
plt.xlabel('$P_{\mathrm{T}}$ / ($\mathrm{bar}$)')

plt.savefig("pareto.png", dpi=600, bbox_inches='tight')
plt.savefig("pareto.pdf",  bbox_inches='tight')
