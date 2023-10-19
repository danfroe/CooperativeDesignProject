# -*- coding: utf-8 -*-
"""
Created on Wed May 31 10:01:59 2023

@author: danfr
"""
import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve
from BenzenePlant_functions import *
import BenzenePlantWithRecycle_StreamData as sd


# PROPERTY DATA:
# Air
molefraction_O2inAir = 0.21 # mol/mol
molefraction_N2inAir = 0.79 # mol/mol
MW_Air = molefraction_N2inAir * MW_N2 + molefraction_O2inAir * MW_O2
weightfraction_O2inAir = molefraction_O2inAir * MW_O2 / MW_Air
# heat capacities + enthalpies at ambient pressure
cp_H2OLiquid25C = 4.182 # kJ/(kg*K)
cp_H2OLiquid100C = 4.216 # kJ/(kg*K)
cp_H2OLiquid_mean = (cp_H2OLiquid25C+cp_H2OLiquid100C)/2
dhv_H2O_100C = 2675.6 - 419.1 #kJ/kg
Tboil_water = 100+273.15 # K
T_ambient = 25 +273.15 # K

#------------------------------------------------------------------------------
# ENERGY BALANCE FURNACE
mHydrogen_stream13 = sd.nHydrogen_stream13 * MW_H2 # kg/h
mMethane_stream13 = sd.nMethane_stream13 * MW_Methane # kg/h
T_fuel = sd.T_stream13                  

Q_H101 = quad(integrand_heatEnthalpyProcessFluid_vapor,sd.T_stream3,sd.T_stream4, \
              args=(sd.massflow_stream3,sd.nHydrogen_stream3, sd.nMethane_stream3, \
                    sd.nBenzene_stream3, sd.nToluene_stream3))

sol = fsolve(energybalancesfurnace, [100,100,100,100], \
             args=(T_ambient,T_fuel, sd.T_fluegasOut,Q_H101,mHydrogen_stream13,mMethane_stream13))
mHydrogen_furnaceIn = sol[0]
mMethane_furnaceIn = sol[1]
mCO2stoichiometric_furnaceOut = sol[2]
mH2Ostoichiometric_furnaceOut = sol[3]
mO2stoichiometric_furnaceIn = 2 * MW_O2/MW_Methane * mMethane_furnaceIn \
                                + 1/2 * MW_O2/MW_H2 * mHydrogen_furnaceIn
mAirstoichiometric_furnaceIn = mO2stoichiometric_furnaceIn / weightfraction_O2inAir
m_fuelgas_furnaceIn = mHydrogen_furnaceIn + mMethane_furnaceIn
print('Heat flow of feed heater H-101: %.2f [kJ/h]' % Q_H101[0])
print('Mass flow of hydrogen into feed heater H-101: %.2f [kg/h]' % mHydrogen_furnaceIn)
print('Mass flow of methane into feed heater H-101: %.2f [kg/h]' % mMethane_furnaceIn)
print('Mass flow of fuel gas into feed heater H-101: %.2f [kg/h]' % m_fuelgas_furnaceIn)
print('Mass flow of oxygen into feed heater H-101: %.2f [kg/h]' % mO2stoichiometric_furnaceIn)
print('Mass flow of air into feed heater H-101: %.2f [kg/h]' % mAirstoichiometric_furnaceIn)
print('Mass flow of CO2 out of feed heater H-101: %.2f [kg/h]' % mCO2stoichiometric_furnaceOut)
print('Mass flow of water out of feed heater H-101: %.2f [kg/h]' % mH2Ostoichiometric_furnaceOut)
