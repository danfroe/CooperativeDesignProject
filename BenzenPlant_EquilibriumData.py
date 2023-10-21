# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 12:23:09 2023

@author: ge37xik
"""

# importing package
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

# create data
t = np.linspace(700,1000,100)

# Define Equilibrium constants from Turton
lnKeq1 = 13.51 + 5037/t - 2.037*np.log(t) + 3.499*10**(-4)*t + 4.173*10**(-8)*t**2 + 3017/t**2 # main reaction
lnKeq2 = 1.788 - 4135.2/t # Side reaction 
y = 0.052+ 0* t


# plot lines
plt.plot(1/t, lnKeq1,  c = 'royalblue', label= "$K_{eq,1}$",linestyle="-",linewidth=2) 
plt.plot(1/t, lnKeq2,  c = 'crimson', label ="$K_{eq,2}$", linestyle="--",linewidth=2)
plt.axvline(x = 1/973, color = 'black', label = '973 K',  linestyle="--",linewidth=2)

# Label axis 
plt.xlabel('$T^{-1}$ /(${\mathrm{K^{-1}}}$)')
plt.ylabel('ln($K_{eq,j}$)/ (-)')

# Layout 
matplotlib.rcParams['axes.linewidth'] = 2
plt.rcParams["font.size"] = "30"
plt.rcParams["figure.figsize"] = (7,5)
plt.rcParams.update({'font.size': 30})

# Plot legend 
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.ylim(-5,10)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

# Save plots 
plt.savefig("eq.png", dpi=600, bbox_inches='tight')
plt.savefig("eq.pdf",  bbox_inches='tight')


t = 973
lnKeq1 = 13.51 + 5037/t - 2.037*np.log(t) + 3.499*10**(-4)*t + 4.173*10**(-8)*t**2 + 3017/t**2
lnKeq2 = 1.788 - 4135.2/t

print(t, np.exp(lnKeq1))
print(t, np.exp(lnKeq2))

plt.show()
