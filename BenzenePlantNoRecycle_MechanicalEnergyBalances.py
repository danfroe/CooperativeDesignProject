# -*- coding: utf-8 -*-
"""
Created on Mon May 22 14:24:32 2023

@author: danfr
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve
from BenzenePlant_functions import *
import BenzenePlantNoRecycle_StreamData as sd

# PROPERTY DATA:
# Densities in kg/kmol
densityToluene50C = 839 # kg/m^3
densityToluene100C = 789.7 # kg/m^3
# Temperatures
T1_lowerToluene = 50 +273.15 # K
T2_upperToluene = 100 +273.15 # K

#------------------------------------------------------------------------------
# ENERGY BALANCES
# P-101: Toluene feed pump
massflow_ToluenePump = sd.massflow_stream1 + sd.massflow_stream9

Cp_Toluene_liquid_stream1 = heatcapacity_interpolation(sd.T_stream1, 50+273.15, 100+273.15, Cp_Toluene_liquid50C, Cp_Toluene_liquid100C)
Cp_Toluene_liquid_stream9 = heatcapacity_interpolation(sd.T_stream9, 50+273.15, 100+273.15, Cp_Toluene_liquid50C, Cp_Toluene_liquid100C)
Treal = fsolve(solveTemperatureOfMixedStreams_liquidBenzeneToluene, 300, \
                          args=(np.array([sd.T_stream1,sd.T_stream9]), 0, np.array([sd.nToluene_stream1,sd.nToluene_stream9]), \
                                0, np.array([Cp_Toluene_liquid_stream1, Cp_Toluene_liquid_stream9])))

densityTolueneReal = density_interpolation(Treal, T1_lowerToluene, T2_upperToluene, densityToluene50C, densityToluene100C)
deltaPressure = sd.p_afterToluenePump - sd.p_beforeToluenePump
w_P101 = ShaftworkOfPump(densityTolueneReal, deltaPressure, 0) # J/kg
W_P101 = w_P101 * massflow_ToluenePump * 1/3600*10**-3 /sd.efficiency_ToluenePump
print('Density of toluene at %i °C: %.2f [kg/m^3]' % (Treal-273.15, densityTolueneReal))
print('Work of toluene feed pump P-101: %.2f [kW]' % W_P101)


# P-102: Reflux pump
massflow_RefluxPump = sd.massflow_condenser
w_P102 = ShaftworkOfPump(1, 0, sd.deltaHeight_RefluxPump) # J/kg; density unimportant but must not be 0
W_P102 = w_P102 * massflow_RefluxPump * 1/3600*10**-3 /sd.efficiency_RefluxPump
print('Work of reflux pump P-102: %.2f [kW]' % W_P102)