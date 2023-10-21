import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.integrate import odeint
from scipy.optimize import fsolve


# initial values for the molar flow rates 1123 K
F_Toluene_0 =  85.5   # kmol/h
F_Hydrogen_0 = 85.5*2  # kmol/h
F_Benzene_0  =   0   # kmol/h
F_Methane_0  =  49.1   # kmol/h

# cold_shots
F_toluene_cold_shot = 0 #kmol/h
F_hydrogen_cold_shot = 29/1.5 #kmol/h
F_benzene_cold_shot = 0 #kmol/h
F_methane_cold_shot = 10 #kmol/h

# Heat capacities
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


# position of the cold shots /m^3 
loc_cs_1 = 0.9
loc_cs_2 = 1.8

# final reactor volume /m^3
loc_end = 3.7


# initial value for the temperature of the reaction mixture 
Trm_0 = 973 # K 
Tcs = 335 #K

# total pressure of the reaction mixture inside the reactor
P_total = 25 # atm 

# initial conditions for all 6 differential equations 
inivalues1=[F_Toluene_0,F_Hydrogen_0,F_Benzene_0,F_Methane_0,Trm_0]

# a function defining our system of differential equations 
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
    
    
    #kinetic parameter for the reaction rate equation from Tsuchiya et al. (1959)
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

# independent variable (volume) array: 150 pts from V=0 to V=60 m^3
Vspan1 = np.linspace(0,loc_cs_1,100) 

# solver output has format [V,[FH,FM,FB,FT,Ts,Trm]]
solver1 = odeint(pfr,inivalues1,Vspan1)


Trm2 = solver1[-1,4]

def alg(z):
    cp_Toluene4,cp_Hydrogen4,cp_Benzene4,cp_Methane4,Trm_new = z
    a = cp_Toluene4 -  cp_Toluene(Trm_new)
    b = cp_Hydrogen4 -  cp_Hydrogen(Trm_new)
    f = cp_Benzene4 -  cp_Benzene(Trm_new)
    d = cp_Methane4 -  cp_Methane(Trm_new)
    e = Trm_new - (solver1[-1,0]*cp_Toluene(Trm2)*Trm2 + solver1[-1,1]*cp_Hydrogen(Trm2)*Trm2 + solver1[-1,2]*cp_Benzene(Trm2)*Trm2 + solver1[-1,3]*cp_Methane(Trm2)*Trm2 + F_toluene_cold_shot*cp_Toluene(Tcs)*Tcs + F_hydrogen_cold_shot*cp_Hydrogen(Tcs)*Tcs + F_benzene_cold_shot*cp_Benzene(Tcs)*Tcs + F_methane_cold_shot*cp_Methane(Tcs)*Tcs)/((F_toluene_cold_shot+solver1[-1,0])*cp_Toluene4 + (F_hydrogen_cold_shot+solver1[-1,1])*cp_Hydrogen4 + (F_benzene_cold_shot+solver1[-1,2])*cp_Benzene4 + (F_methane_cold_shot+solver1[-1,3])*cp_Methane4)
    return(a,b,f,d,e)    

solution = fsolve(alg, [3,3,3,3,3])

Trm_new = solution[4]
#print("Temperature nach cold_shot_1",Trm_new)

inivalues2 = [solver1[-1,0]+F_toluene_cold_shot,solver1[-1,1]+F_hydrogen_cold_shot,solver1[-1,2]+F_benzene_cold_shot,solver1[-1,3]+F_methane_cold_shot,Trm_new]

# independent variable (volume) array: 150 pts from V=0 to V=60 m^3
Vspan2 = np.linspace(loc_cs_1,loc_cs_2,100) 

# solver output has format [V,[FH,FM,FB,FT,Ts,Trm]]
solver2 = odeint(pfr,inivalues2,Vspan2)

Trm3 = solver2[-1,4]


def alg2(z):
    cp_Toluene5,cp_Hydrogen5,cp_Benzene5,cp_Methane5,Trm_new2 = z
    a = cp_Toluene5 -  cp_Toluene(Trm_new2)
    b = cp_Hydrogen5 -  cp_Hydrogen(Trm_new2)
    f = cp_Benzene5 -  cp_Benzene(Trm_new2)
    d = cp_Methane5 -  cp_Methane(Trm_new2)
    e = Trm_new2 - (solver2[-1,0]*cp_Toluene(Trm3)*Trm3 + solver2[-1,1]*cp_Hydrogen(Trm3)*Trm3 + solver2[-1,2]*cp_Benzene(Trm3)*Trm3 + solver2[-1,3]*cp_Methane(Trm3)*Trm3 + F_toluene_cold_shot*cp_Toluene(Tcs)*Tcs + F_hydrogen_cold_shot*cp_Hydrogen(Tcs)*Tcs + F_benzene_cold_shot*cp_Benzene(Tcs)*Tcs + F_methane_cold_shot*cp_Methane(Tcs)*Tcs)/((F_toluene_cold_shot+solver2[-1,0])*cp_Toluene5 + (F_hydrogen_cold_shot+solver2[-1,1])*cp_Hydrogen5 + (F_benzene_cold_shot+solver2[-1,2])*cp_Benzene5 + (F_methane_cold_shot+solver2[-1,3])*cp_Methane5)
    return(a,b,f,d,e)    

solution2 = fsolve(alg2, [300,300,300,300,1113])

Trm_new2 = solution2[4]


inivalues3 = [solver2[-1,0]+F_toluene_cold_shot,solver2[-1,1]+F_hydrogen_cold_shot,solver2[-1,2]+F_benzene_cold_shot,solver2[-1,3]+F_methane_cold_shot,Trm_new2]

# independent variable (volume) array: 150 pts from V=0 to V=60 m^3
Vspan3 = np.linspace(loc_cs_2,loc_end,100) 

# solver output has format [V,[FH,FM,FB,FT,Ts,Trm]]
solver3 = odeint(pfr,inivalues3,Vspan3)


# concatenate the results of solver1, solver2 and solver3 
toluene=np.concatenate((solver1[:,0], solver2[:,0],solver3[:,0]))
hydrogen= np.concatenate((solver1[:,1], solver2[:,1],solver3[:,1]))
benzene=np.concatenate((solver1[:,2], solver2[:,2],solver3[:,2]))
methane=np.concatenate((solver1[:,3], solver2[:,3],solver3[:,3]))
temperature=np.concatenate((solver1[:,4], solver2[:,4],solver3[:,4]))

Span=np.concatenate((Vspan1,Vspan2,Vspan3))

#Figure options 
fig, (ax1, ax2) = plt.subplots(2)
fig.tight_layout(pad=0.3)
fig = matplotlib.pyplot.gcf()
fig.set_size_inches(7, 7)
plt.rcParams["font.size"] = "24"
plt.setp(ax1.spines.values(), linewidth=1.7)
plt.setp(ax2.spines.values(), linewidth=1.7)

# Plot commands
ax1.plot(Span,toluene, label="T" , c="r")
ax1.plot(Span,hydrogen, label="H", c="c", linestyle="-")
ax1.plot(Span,methane,label="M", c="blue", linestyle="-")
ax1.plot(Span,benzene,label="B", c="forestgreen", linestyle="-")
ax2.plot(Span,temperature, c ='red')



# axis commands
ax1.set(xlabel = '$V_{\mathrm{R}}$ / (${\mathrm{m^3}}$)')
ax1.set(ylabel = '$F^{\mathrm{i}}$ / (kmol/h)')
ax1.set_xlim([0,3.6])
ax1.set_ylim([0, 200])
ax2.set_xlim([0,3.6])
ax1.set_xticks([0,  0.9,  1.8,  2.7, 3.6])
ax2.set_xticks([0,  0.9,  1.8,  2.7, 3.6])
ax2.set(xlabel = '$V_{\mathrm{R}}$ / (${\mathrm{m^3}}$)')
ax2.set(ylabel = '$T_{\mathrm{m}}$ / (K)')
ax1.xaxis.set_tick_params(width=2)
ax1.yaxis.set_tick_params(width=2)
ax2.xaxis.set_tick_params(width=2)
ax2.yaxis.set_tick_params(width=2)

# Legend
ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))

# Temperature range
ax2.axhline(y=1025, color='gray',linestyle="--")
ax2.axhline(y=973, color='gray',linestyle="--")
ax2.set_ylim([950,1050])

# Calculation of conversion as a function of volume
X_target=solver1[:,3]
X_toluene1=[]
for i in range(100):
    X_toluene1.append((F_Toluene_0-X_target[i])/F_Toluene_0) #(intial-final)/initial
    
X_target=solver3[:,0]
X=[]
for i in range(100):
    X.append((F_Toluene_0-X_target[i])/F_Toluene_0) #(intial-final)/initial
volumeindex = 99
for i in range(100):
    if X[i] >= 0.75 and i <= volumeindex:
        volumeindex = i 
        
# Print converion and total reactor volume         
#print("Total reactor volume:", Vspan3[volumeindex])
#print(solver3[volumeindex,0])

# Save Figures
plt.savefig("simulation file 1.png", dpi=600, bbox_inches='tight')
plt.savefig("simulation file 1.pdf",  bbox_inches='tight')

# Show figure
plt.show()
