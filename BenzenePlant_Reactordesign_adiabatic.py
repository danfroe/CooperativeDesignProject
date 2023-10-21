import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.integrate import odeint


# Define cp-values as a function of temperature 

def cp_Toluene(Trm): #kJ/kmol/K
    cp_Toluene= (-2.435*10 + 5.125*0.1*Trm + -2.765e-4* Trm**2 +4.911e-8*Trm**3)
    return cp_Toluene

def cp_Hydrogen(Trm):#kJ/kmol/K
    cp_Hydrogen =  (3.224 * 10 + 1.924 * 10**-3*Trm + 1.055 * 10**-5* Trm**2 + -3.596*10**-9*Trm**3)
    return cp_Hydrogen

def cp_Benzene(Trm):#kJ/kmol/K
    cp_Benzene=( -3.392 * 10+ 4.739* 0.1*Trm + -3.017 * 10**-4* Trm**2 +7.130*10**-8 *Trm**3)
    return cp_Benzene

def cp_Methane(Trm):#kJ/kmol/K
    cp_Methane=( 1.925 * 10+ 5.213* 10**-2*Trm + 1.197 * 10**-5* Trm**2 + -1.132*10**-8*Trm**3)
    return cp_Methane



# Initial values for the molar flow rates 
F_Toluene_0 =  85.5   # kmol/h
F_Hydrogen_0 = 85.5*2   # kmol/h
F_Benzene_0  =   0   # kmol/h
F_Methane_0  =  49.1   # kmol/h



# Initial value for the temperature of the reaction mixture 
Trm_0 = 973 # K 


# total pressure of the reaction mixture inside the reactor
P_total = 25 # atm 



# Initial conditions for all 6 differential equations 
inivalues=[F_Toluene_0,F_Hydrogen_0,F_Benzene_0,F_Methane_0,Trm_0]


# A function defining our system of differential equations 
def pfr(z,V):
    
    [F_Toluene,F_Hydrogen,F_Benzene,F_Methane,Trm] = z
    
    F_total = F_Hydrogen + F_Methane + F_Benzene + F_Toluene #Total molar flowrate

    # enthalpy of reaction as a function of temperature kj/kmol 
    DHrxn = - 37190 - 17.24 * Trm + 29.09 * 10 ** (-4) * Trm**2 + 0.6939 * 10**(-6)*Trm**3+50160/Trm
    
    x_Hydrogen = F_Hydrogen / F_total  # mole fractions 
    x_Methane  = F_Methane  / F_total  # mole fractions 
    x_Benzene  = F_Benzene  / F_total  # mole fractions 
    x_Toluene  = F_Toluene  / F_total  # mole fractions 
    
    p_Hydrogen = x_Hydrogen * P_total  # partial pressures 
    p_Methane  = x_Methane  * P_total  # partial pressures
    p_Benzene  = x_Benzene  * P_total  # partial pressures
    p_Toluene  = x_Toluene  * P_total  # partial pressures
    
    # kinetic parameter for the reaction rate equation 
    k = 3.55 * 10**5 * np.exp(-54500/1.987/Trm) * 3600 * 1000 # kmol/m^3/h/atm^1.5
    
    # definition of the reaction rate equation
    r_Benzene =  k * p_Hydrogen**0.5 * p_Toluene  # kmol/m^3/h
    
    # definition of the system of differential equations 
    dF_HydrogendV  = - r_Benzene
    dF_MethanedV   = + r_Benzene
    dF_BenzenedV   = + r_Benzene
    dF_ToluenedV   = - r_Benzene
    dTrm_dV = ( -DHrxn * r_Benzene ) / (F_Toluene*cp_Toluene(Trm)+F_Hydrogen*cp_Hydrogen(Trm)+F_Benzene*cp_Benzene(Trm)+F_Methane*cp_Methane(Trm))

    
    return dF_ToluenedV, dF_HydrogendV, dF_BenzenedV, dF_MethanedV, dTrm_dV

# Independent variable (volume) array: 150 pts from V=0 to V=60 m^3
Vspan = np.linspace(0,2.5,150) 

# Solver output has format [V,[FH,FM,FB,FT,Ts,Trm]]
solver = odeint(pfr,inivalues,Vspan)


# Calculation of conversion as a function of volume
X_target=solver[:,0]
X=[]
for i in range(150):
    X.append((F_Toluene_0-X_target[i])/F_Toluene_0) #(intial-final)/initial
    
volumeindex = 149
for i in range(150):
    if X[i] >= 0.75 and i <= volumeindex:
        volumeindex = i 
print("Reactor volume:", Vspan[volumeindex])

print("Final temperature:", solver[-1,4])

# Figure options 
fig, (ax1, ax2)= plt.subplots(2)
fig.tight_layout(pad=0.3)
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(7, 7)
plt.rcParams.update({'font.size': 20})

# Plot commands
ax1.plot(Vspan,solver[:,0], label='T', c="red")
ax1.plot(Vspan,solver[:,1], label='H', c="c")
ax1.plot(Vspan,solver[:,2], label='B', c="forestgreen")
ax1.plot(Vspan,solver[:,3], label='M', c = 'blue')

ax2.plot(Vspan,solver[:,4], c ='red')


# Plot Conversion 
#ax3.plot(Vspan,X, label= '$X_{Toluene}$', c="magenta")
#ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))
#ax3.set_yticks([0,0.5,1])
#ax3.set(ylabel = '$X_{Toluene}$ /(mol/mol)')
#ax3.axhline(y=0.75, color='black',linestyle="--")


# axis commands
ax1.set(xlabel = '$V$ /($m^3$)')
ax1.set(ylabel = '$F_i$ /(kmol/h)')
ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax1.set_ylim([0, 200])

ax2.set(xlabel = '$V$ /($m^3$)')
ax2.set(ylabel = '$T_j$ /(K)')
ax2.set_ylim([950,1125])



#plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
# show figure
plt.show()