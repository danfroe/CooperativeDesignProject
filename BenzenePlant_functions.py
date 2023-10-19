# -*- coding: utf-8 -*-
"""
Created on Thu May 25 10:05:26 2023

@author: danfr
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve


# PROPERTY DATA
R = 8.314 # kJ/(kmol*K), gas constant
g = 9.81 # m/s^2
# Molar weights in kg/kmol
MW_H2 = 2.016
MW_H2O = 18.015
MW_Methane = 16.0436
MW_Benzene = 78.113
MW_Toluene = 92.14
MW_O2 = 32
MW_CO2 = 44
MW_N2 = 28
# heat capacities
Cp_Benzene_liquid50C = 1.802 # kJ/(kg*K), from VDI Heat Atlas
Cp_Benzene_liquid100C = 1.980 # kJ/(kg*K)
Cp_Benzene_liquid150C = 2.189 # kJ/(kg*K)
Cp_Toluene_liquid50C = 1.783 # kJ/(kg*K)
Cp_Toluene_liquid100C = 1.964 # kJ/(kg*K)
Cp_Toluene_liquid150C = 2.158 # kJ/(kg*K)
Cp_Toluene_liquid200C = 2.380 # kJ/(kg*K)
Cp_Toluene_liquid250C = 2.697 # kJ/(kg*K)
# Lower heating values
LHV_Hydrogen = 120000 # kJ/kg, from engineeringtoolbox.com
LHV_Methane = 50000 # kJ/kg, from engineeringtoolbox.com

#------------------------------------------------------------------------------
# FUNCTIONS:
# Heat capacities in kJ/(kg*K)
def Cp_Steam(T):
    return (3.224e1 + 1.924e-3*T + 1.055e-5*T**2 + (-3.596e-9)*T**3)/MW_H2O

def Cp_Hydrogen(T):
    return (3.224e1 + 1.924e-3*T + 1.055e-5*T**2 + (-3.596e-9)*T**3)/MW_H2

def Cp_Methane(T):
    return (1.925e1 + 5.213e-2*T + 1.197e-5*T**2 + (-1.132e-8)*T**3)/MW_Methane

def Cp_Benzene_vapor(T):
    return (-3.392e1 + 4.739e-1*T + (-3.017e-4)*T**2 + (7.130e-8)*T**3)/MW_Benzene

def Cp_Toluene_vapor(T):
    return (-2.435e1 + 5.125e-1*T + (-2.765e-4)*T**2 + 4.911e-8*T**3)/MW_Toluene

def Cp_CO2(T):
    return (1.98e1 + 7.344e-2*T + (-5.602e-5)*T**2 + 1.715e-8*T**3)/MW_CO2

def Cp_processfluid_vapor(T,nHydrogen,nMethane,nBenzene,nToluene): # insert values Hy,Me,Be,To in kmol/h
    TotalMassFlow = MassflowProcessfluid(nHydrogen,nMethane,nBenzene,nToluene)
    xHydrogen = nHydrogen*MW_H2/TotalMassFlow
    xMethane = nMethane*MW_Methane/TotalMassFlow
    xBenzene = nBenzene*MW_Benzene/TotalMassFlow
    xToluene = nToluene*MW_Toluene/TotalMassFlow
    return xHydrogen*Cp_Hydrogen(T) + xMethane*Cp_Methane(T) + xBenzene*Cp_Benzene_vapor(T) + xToluene*Cp_Toluene_vapor(T)

def Cp_processfluid_liquidBenzeneToluene(nBenzene, nToluene, Cp_Benzene_liquid, Cp_Toluene_liquid): # insert values Hy,Me,Be,To in kmol/h
    TotalMassFlow = MassflowProcessfluid(0, 0, nBenzene, nToluene)
    xBenzene = nBenzene*MW_Benzene/TotalMassFlow
    xToluene = nToluene*MW_Toluene/TotalMassFlow
    return xBenzene*Cp_Benzene_liquid + xToluene*Cp_Toluene_liquid

def MassflowProcessfluid(nHydrogen,nMethane,nBenzene,nToluene):
    massflow = nHydrogen*MW_H2 + nMethane*MW_Methane + nBenzene*MW_Benzene + nToluene*MW_Toluene
    return massflow

def enthalpy_interpolation(Treal,Tlower,Tupper,Hlower,Hupper): # T in °C or K
    return Hlower + (Hupper-Hlower)/(Tupper-Tlower) * (Treal-Tlower)

def heatcapacity_interpolation(Treal,Tlower,Tupper,Cplower,Cpupper): # T in °C or K
    return Cplower + (Cpupper-Cplower)/(Tupper-Tlower) * (Treal-Tlower)

def dhv_mixture(dhv_Benzene,dhv_Toluene,nHydrogen,nMethane,nBenzene,nToluene):
    TotalMassFlow = MassflowProcessfluid(nHydrogen,nMethane,nBenzene,nToluene)
    xBenzene = nBenzene*MW_Benzene/TotalMassFlow
    xToluene = nToluene*MW_Toluene/TotalMassFlow
    return xBenzene*dhv_Benzene + xToluene*dhv_Toluene

def dhv_WatsonsEquation(T1, T2, T_crit, dhv1): # from Himmelblau
    Tr1 = T1 / T_crit
    Tr2 = T2 / T_crit
    dhv2 = ((1-Tr2)/(1-Tr1))**0.38 * dhv1
    return dhv2

def solveTemperatureOfMixedStreams_vapor(T, array_T_in, array_nHydrogen_in, array_nMethane_in, array_nBenzene_in, array_nToluene_in): # input as 0 or np.array
    dh_streamsIn = (array_nHydrogen_in*MW_H2*Cp_Hydrogen(array_T_in) + array_nMethane_in*MW_Methane*Cp_Methane(array_T_in) \
        + array_nBenzene_in*MW_Benzene*Cp_Benzene_vapor(array_T_in) + array_nToluene_in*MW_Toluene*Cp_Toluene_vapor(array_T_in))*array_T_in
    dh_streamUnknown = (np.sum(array_nHydrogen_in)*MW_H2*Cp_Hydrogen(T) + np.sum(array_nMethane_in)*MW_Methane*Cp_Methane(T) \
           + np.sum(array_nBenzene_in)*MW_Benzene*Cp_Benzene_vapor(T) + np.sum(array_nToluene_in)*MW_Toluene*Cp_Toluene_vapor(T))*T
    return np.sum(dh_streamsIn) - dh_streamUnknown

def solveTemperatureOfMixedStreams_liquidBenzeneToluene(T, array_T_in, array_nBenzene_in, array_nToluene_in, array_Cp_Benzene_liquid, array_Cp_Toluene_liquid): # input as 0 or np.array
    dh_streamsIn = (array_nBenzene_in*MW_Benzene*array_Cp_Benzene_liquid + array_nToluene_in*MW_Toluene*array_Cp_Toluene_liquid)*array_T_in
    dh_streamUnknown = (+ np.sum(array_nBenzene_in)*MW_Benzene*heatcapacity_interpolation(T, 50+273.15, 100+273.15, Cp_Benzene_liquid50C, Cp_Benzene_liquid100C) \
                        + np.sum(array_nToluene_in)*MW_Toluene*heatcapacity_interpolation(T, 50+273.15, 100+273.15, Cp_Toluene_liquid50C, Cp_Toluene_liquid100C))*T
    return np.sum(dh_streamsIn) - dh_streamUnknown

def integrand_heatEnthalpyProcessFluid_vapor(T,massflow_processfluid,nHydrogen,nMethane,nBenzene,nToluene):
    cp = Cp_processfluid_vapor(T,nHydrogen,nMethane,nBenzene,nToluene) # kJ/(kg*K)
    dH = massflow_processfluid*cp
    return dH

def massflow_coolingwater(dH,cp_water,T_cw_in,T_cw_out):
    return - dH / (cp_water*(T_cw_out-T_cw_in))

def density_interpolation(Treal,T1_lower,T2_upper,densityT1,densityT2): # T in °C or K
    return densityT1 + (densityT2-densityT1)/(T2_upper-T1_lower) * (Treal-T1_lower)
    
def ShaftworkOfPump(density,DeltaPressure,DeltaHeight):
    Shaftwork = DeltaPressure / density + g * DeltaHeight # friction neglected, no velocity change; from Introduction to Process Engineering script
    return Shaftwork

def ShaftworkOfCompressor(var, Cp, T_in, p_in, p_out):
    Shaftwork = var[0]
    T_out = var[1]
    return [Shaftwork - Cp * (T_out - T_in),
            Shaftwork - Cp * T_in * ((p_out/p_in)**(R/Cp) - 1)] # from "Perry's Chemical Engineer Handbook" p. 4-7

def Antoine_equation_Toluene(P): # P in bar; from NIST/Ambrose
    A = 4.54436
    B = 1738.123
    C = 0.394
    T = -B / (np.log10(P) - A) - C
    return T


# HEATER
def integrand_heatEnthalpyFuel(T,mHydrogen,mMethane):
    dh = mHydrogen * Cp_Hydrogen(T) + mMethane * Cp_Methane(T)
    return dh

def integrand_heatEnthalpyFlueGasout(T,mCO2,mH2O):
    dh = mCO2 * Cp_CO2(T) + mH2O * Cp_Steam(T)
    return dh

        
def energybalancesfurnace(var,T_ambient,T_fuel,T_fluegasOut,Q_H101,mHydrogen_available,mMethane_available):
    mHydrogen = var[0]
    mMethane = var[1]
    mCO2 = var[2]
    mH2O = var[3]
    dh_heatFuel = quad(integrand_heatEnthalpyFuel,T_ambient,T_fuel,args=(mHydrogen,mMethane))
    dh_heatFlueGasout = quad(integrand_heatEnthalpyFlueGasout,T_ambient,T_fluegasOut,args=(mCO2,mH2O))
    return [
            mHydrogen*LHV_Hydrogen + mMethane*LHV_Methane + dh_heatFuel[0] - dh_heatFlueGasout[0] - Q_H101[0],
            mMethane_available/mHydrogen_available * mHydrogen - mMethane,
            MW_CO2/MW_Methane * mMethane - mCO2,
            2 * MW_H2O/MW_Methane * mMethane + MW_H2O/MW_H2 * mHydrogen - mH2O
            ]