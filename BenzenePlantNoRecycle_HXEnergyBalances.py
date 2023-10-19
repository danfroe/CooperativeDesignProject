# -*- coding: utf-8 -*-
"""
Created on Thu May  4 16:39:24 2023

@author: danfr
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve
from BenzenePlant_functions import *
import BenzenePlantNoRecycle_StreamData as sd

# PROPERTY DATA:
cp_water = 4.18 # cp of liquid water in kJ/(kg*K)
# Enthalpies of vaporization/condensation
dhv_Benzene100C = 378.6 # kJ/kg
dhv_Benzene150C = 335.9 # kJ/kg
dhv_Toluene100C = 367.7 # kJ/kg
dhv_Toluene150C = 333.0 # kJ/kg
# T critical
T_crit_Benzene = 562 # K, from Engineering Toolbox
T_crit_Toluene = 592 # K, from Leandros' greek book / Engineering Toolbox

#------------------------------------------------------------------------------
# UTILITIES:
# Low pressure steam (5 barg)
dhc_lps = 675.57 - 2757.4 # kJ/kg; from VDI Heat Atlas p. 178 (160 °C, 6.18 bar)
print('Enthalpy of condensation of lps: %.2f [kJ/kg]' % dhc_lps)
# Medium pressure steam (10 barg)
dhc_mps = enthalpy_interpolation(sd.T_mps,180+273.15,190+273.15,(763.19-2777.2),(807.57-2785.3)) # kJ/kg; from VDI Heat Atlas p. 178
print('Enthalpy of condensation of mps: %.2f [kJ/kg]' % dhc_mps)
# High pressure steam (41 barg)
dhc_hps = enthalpy_interpolation(sd.T_hps,250+273.15,260+273.15,(1085.7-2801.0),(1134.8-2796.6)) # kJ/kg; from VDI Heat Atlas p. 178
print('Enthalpy of condensation of hps: %.2f [kJ/kg]' % dhc_hps)
print()


# ENERGY BALANCES:
# Heat exchanger E-101: Feed preheater
# streams 1 + 9 + 11 + 12 --> 2
Tboil_Toluene = sd.T_stream3
T_stream2vapors = fsolve(solveTemperatureOfMixedStreams_vapor, 300, \
                         args=(sd.T_stream11, sd.nHydrogen_stream11-sd.nHydrogen_stream10, 0, 0, 0))

Cp_Toluene_liquid_stream1 = heatcapacity_interpolation(sd.T_stream1, 50+273.15, 100+273.15, Cp_Toluene_liquid50C, Cp_Toluene_liquid100C)
Cp_Toluene_liquid_stream9 = heatcapacity_interpolation(sd.T_stream9, 50+273.15, 100+273.15, Cp_Toluene_liquid50C, Cp_Toluene_liquid100C)
T_stream2liquids = fsolve(solveTemperatureOfMixedStreams_liquidBenzeneToluene, 300, \
                          args=(np.array([sd.T_stream1,sd.T_stream9]), 0, np.array([sd.nToluene_stream1,sd.nToluene_stream9]), \
                                0, np.array([Cp_Toluene_liquid_stream1, Cp_Toluene_liquid_stream9])))
#---------
nHydrogen_in = np.array([0, 0, sd.nHydrogen_stream11-sd.nHydrogen_stream10])
nToluene_in = np.array([sd.nToluene_stream1, sd.nToluene_stream9, 0])
T_in = np.array([sd.T_stream1, sd.T_stream9, sd.T_stream11])
T_stream2 = fsolve(solveTemperatureOfMixedStreams_vapor, 300, args=(T_in,nHydrogen_in,0,0,nToluene_in))
#---------
massflow_stream2vapors = MassflowProcessfluid(sd.nHydrogen_stream11-sd.nHydrogen_stream10, 0, 0, 0)
massflow_stream2liquids = MassflowProcessfluid(0, 0, 0, sd.nToluene_stream1+sd.nToluene_stream9)

Q_E101vapors = quad(integrand_heatEnthalpyProcessFluid_vapor,T_stream2vapors[0],sd.T_stream3, \
               args=(massflow_stream2vapors, sd.nHydrogen_stream11-sd.nHydrogen_stream10, 0, 0, 0))

Cp_Tolueneliquid = heatcapacity_interpolation((T_stream2liquids[0]+Tboil_Toluene)/2, 200+273.15, 250+273.15, Cp_Toluene_liquid200C, Cp_Toluene_liquid250C)
dhv_Toluene = dhv_WatsonsEquation(150+273.15, Tboil_Toluene, T_crit_Toluene, dhv_Toluene150C)
Q_E101liquids = massflow_stream2liquids * Cp_Tolueneliquid * (Tboil_Toluene - T_stream2liquids[0]) \
                + massflow_stream2liquids * dhv_Toluene

Q_E101 = Q_E101vapors + Q_E101liquids
massflow_hps_E101 = -Q_E101[0]/dhc_hps
print('T_in of vapors of HX E-101: %.2f [K]' % T_stream2vapors[0])
print('T_in of liquids of HX E-101: %.2f [K]' % T_stream2liquids[0])
print('T_in overall of HX E-101: %.2f [K]' % T_stream2[0])
print('Boiling point of toluene: %.2f [K]' % Tboil_Toluene)
print('Heat flow of feed preheater E-101: %.2f [kJ/h]' % Q_E101[0])
print('Mass flow of high pressure steam in feed preheater E-101: %.2f [kg/h]' % massflow_hps_E101)
print()


# Heat exchanger E-102: Reactor effluent cooler
Q_E102 = quad(integrand_heatEnthalpyProcessFluid_vapor,sd.T_stream5,sd.T_stream6, \
              args=(sd.massflow_stream5,sd.nHydrogen_stream5, sd.nMethane_stream5, \
                    sd.nBenzene_stream5, sd.nToluene_stream5))
massflow_cw_E102 = massflow_coolingwater(Q_E102[0],cp_water,sd.T_cw_in,sd.T_cw_out) # kg/h
print('Heat flow of reactor effluent cooler E-102: %.2f [kJ/h]' % Q_E102[0])
print('Mass flow of cooling water in reactor effluent cooler E-102: %.2f [kg/h]' % massflow_cw_E102)
print()


# Heat exchanger E-103: Tower feed heater
Cp_Benzene_liquid = heatcapacity_interpolation((sd.T_stream6+sd.T_stream7)/2, 50+273.15, 100+273.15, Cp_Benzene_liquid50C, Cp_Benzene_liquid100C)
Cp_Toluene_liquid = heatcapacity_interpolation((sd.T_stream6+sd.T_stream7)/2, 50+273.15, 100+273.15, Cp_Toluene_liquid50C, Cp_Toluene_liquid100C)
Q_E103 = sd.massflow_stream7 * Cp_processfluid_liquidBenzeneToluene(sd.nBenzene_stream7, sd.nToluene_stream7, Cp_Benzene_liquid, Cp_Toluene_liquid) * (sd.T_stream7-sd.T_stream6)
massflow_lps_E103 = -Q_E103/dhc_lps
print('Heat flow of tower feed heater E-103: %.2f [kJ/h]' % Q_E103)
print('Mass flow of low pressure steam in tower feed heater E-103: %.2f [kg/h]' % massflow_lps_E103)
print()


# Heat exchanger E-104: Benzene condenser
Q_E104_vapor = quad(integrand_heatEnthalpyProcessFluid_vapor,sd.T_condenser_in,sd.T_condenser_out, \
              args=(sd.massflow_condenser,sd.nHydrogen_condenser, sd.nMethane_condenser, \
                    sd.nBenzene_condenser, sd.nToluene_condenser))
dhv_Benzene = dhv_WatsonsEquation(150+273.15, sd.T_condenser_out, T_crit_Benzene, dhv_Benzene150C)
dhv_Toluene = dhv_WatsonsEquation(150+273.15, sd.T_condenser_out, T_crit_Toluene, dhv_Toluene150C)
Q_E104_cond = sd.massflow_condenser*dhv_mixture(-dhv_Benzene, -dhv_Toluene, sd.nHydrogen_condenser, \
                                        sd.nMethane_condenser, sd.nBenzene_condenser, sd.nToluene_condenser)
Q_E104 = Q_E104_vapor[0] + Q_E104_cond
massflow_cw_E104 = massflow_coolingwater(Q_E104,cp_water,sd.T_cw_in,sd.T_cw_out)
print('Heat flow of benzene condenser E-104: %.2f [kJ/h]' % Q_E104)
print('Mass flow of cooling water in benzene condenser E-104: %.2f [kg/h]' % massflow_cw_E104)
print()


# Heat exchanger E-105: Product cooler
Cp_Benzene_liquid = heatcapacity_interpolation((sd.T_condenser_out+sd.T_stream8)/2, 50+273.15, 100+273.15, Cp_Benzene_liquid50C, Cp_Benzene_liquid100C)
Q_E105 = sd.massflow_stream8 * Cp_Benzene_liquid * (sd.T_stream8-sd.T_condenser_out)
massflow_cw_E105 = massflow_coolingwater(Q_E105,cp_water,sd.T_cw_in,sd.T_cw_out)
print('Heat flow of product cooler E-105: %.2f [kJ/h]' % Q_E105)
print('Mass flow of cooling water in product cooler E-105: %.2f [kg/h]' % massflow_cw_E105)
print()


# Heat exchanger E-106: Benzene reboiler
#Cp_Benzene_liquid_stream7 = heatcapacity_interpolation(sd.T_stream7, 50+273.15, 100+273.15, Cp_Benzene_liquid50C, Cp_Benzene_liquid100C)
#Cp_Benzene_liquid_recycle = heatcapacity_interpolation(sd.T_condenser, 100+273.15, 150+273.15, Cp_Benzene_liquid100C, Cp_Benzene_liquid150C)
#Cp_Toluene_liquid_stream7 = heatcapacity_interpolation(sd.T_stream7, 50+273.15, 100+273.15, Cp_Toluene_liquid50C, Cp_Toluene_liquid100C)
#T_beforereboiler = fsolve(solveTemperatureOfMixedStreams_liquidBenzeneToluene, 30, \
#                          args=(np.array([sd.T_stream7, sd.T_condenser]), np.array([sd.nBenzene_stream7, sd.nBenzene_recycle]), \
#                                np.array([sd.nToluene_stream7, 0]), np.array([Cp_Benzene_liquid_stream7, Cp_Benzene_liquid_recycle]), \
#                                np.array([Cp_Toluene_liquid_stream7, 0])))

dhv_Benzene = dhv_WatsonsEquation(150+273.15, sd.T_reboiler, T_crit_Benzene, dhv_Benzene150C)
dhv_Toluene = dhv_WatsonsEquation(150+273.15, sd.T_reboiler, T_crit_Toluene, dhv_Toluene150C)
#Cp_Benzene_liquid = heatcapacity_interpolation((T_beforereboiler+sd.T_reboiler)/2, 100+273.15, 150+273.15, Cp_Benzene_liquid100C, Cp_Benzene_liquid150C)
#Cp_Toluene_liquid = heatcapacity_interpolation((T_beforereboiler+sd.T_reboiler)/2, 100+273.15, 150+273.15, Cp_Toluene_liquid100C, Cp_Toluene_liquid150C)

#Q_E106 = sd.massflow_reboiler * Cp_processfluid_liquidBenzeneToluene(sd.nBenzene_reboiler, sd.nToluene_reboiler, Cp_Benzene_liquid, Cp_Toluene_liquid) * (sd.T_reboiler-T_beforereboiler) \
#        + sd.massflow_reboiler*dhv_mixture(dhv_Benzene, dhv_Toluene, 0,  0, sd.nBenzene_reboiler, sd.nToluene_reboiler)
Q_E106 = sd.massflow_reboiler*dhv_mixture(dhv_Benzene, dhv_Toluene, 0,  0, sd.nBenzene_reboiler, sd.nToluene_reboiler)
massflow_mps_E106 = -Q_E106/dhc_mps
print('Heat flow of benzene reboiler E-106: %.2f [kJ/h]' % Q_E106)
print('Mass flow of medium pressure steam in benzene reboiler E-106: %.2f [kg/h]' % massflow_mps_E106)
