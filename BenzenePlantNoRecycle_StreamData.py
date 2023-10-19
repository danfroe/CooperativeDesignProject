# -*- coding: utf-8 -*-
"""
Created on Thu Jun  1 14:02:44 2023

@author: danfr
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve
from BenzenePlant_functions import *

# PUMPS/COMPRESSOR
p_beforeToluenePump = 2*10**5 # Pascal
p_afterToluenePump = 25*10**5 # Pascal
efficiency_ToluenePump = 0.75 # in Turton: 0.75

deltaHeight_RefluxPump = 8+3 # m
efficiency_RefluxPump = 0.5 # in Turton: 0.5


# STREAM DATA
T_stream1 = 25 +273.15 # K
nToluene_stream1 = 64.1 # kmol/h
massflow_stream1 = MassflowProcessfluid(0, 0, 0, nToluene_stream1)

nToluene_stream2 = 85.5 # kmol/h
nHydrogen_stream2 = 170.9 # kmol/h
nBenzene_stream2 = 0 # kmol/h
nMethane_stream2 = 0# kmol/h

# T_stream3 is the boiling temperature of toluene, see HX Energy Balance
T_stream3 = Antoine_equation_Toluene(p_afterToluenePump/10**5)
nToluene_stream3 = nToluene_stream2
nHydrogen_stream3 = nHydrogen_stream2
nBenzene_stream3 = nBenzene_stream2
nMethane_stream3 = nMethane_stream2
massflow_stream3 = MassflowProcessfluid(nHydrogen_stream3, nMethane_stream3, \
                                        nBenzene_stream3, nToluene_stream3)

T_stream4 = 973 # K; from Steffi 973 K
nToluene_stream4 = nToluene_stream2
nHydrogen_stream4 = nHydrogen_stream2
nBenzene_stream4 = nBenzene_stream2
nMethane_stream4 = nMethane_stream2

T_stream5 = 1021 # K; from Steffi
nToluene_stream5 = 21.4 # kmol/h
nHydrogen_stream5 = 106.8 # kmol/h
nBenzene_stream5 = 64.1 # kmol/h
nMethane_stream5 = 64.1 # kmol/h
massflow_stream5 = MassflowProcessfluid(nHydrogen_stream5, nMethane_stream5, \
                                        nBenzene_stream5, nToluene_stream5)

T_stream6 = 38 +273.15 # K; temporary value from Turton
nToluene_stream6 = 0 # kmol/h
nHydrogen_stream6 = nHydrogen_stream5 # kmol/h
nBenzene_stream6 = 0 # kmol/h
nMethane_stream6 = nMethane_stream5 # kmol/h

T_stream7 = 126.65 +273.15 # K; temporary value from Turton
nToluene_stream7 = 21.4 # kmol/h
nHydrogen_stream7 = 0 # kmol/h
nBenzene_stream7 = 64.1 # kmol/h
nMethane_stream7 = 0 # kmol/h
massflow_stream7 = MassflowProcessfluid(nHydrogen_stream7, nMethane_stream7, \
                                        nBenzene_stream7, nToluene_stream7)

T_stream8 = 38 +273.15 # K; temporary value from Turton
nToluene_stream8 = 0 # kmol/h
nHydrogen_stream8 = 0 # kmol/h
nBenzene_stream8 = 64.1 # kmol/h
nMethane_stream8 = 0 # kmol/h
massflow_stream8 = MassflowProcessfluid(nHydrogen_stream8, nMethane_stream8, \
                                        nBenzene_stream8, nToluene_stream8)

T_stream9 = 147 +273.15 # K; temporary value from Turton
nToluene_stream9 = 21.4 # kmol/h
nHydrogen_stream9 = 0 # kmol/h
nBenzene_stream9 = 0 # kmol/h
nMethane_stream9 = 0 # kmol/h
massflow_stream9 = MassflowProcessfluid(nHydrogen_stream9, nMethane_stream9, \
                                        nBenzene_stream9, nToluene_stream9)

T_stream10 = 25 +273.15 # K
nToluene_stream10 = 0 # kmol/h
nHydrogen_stream10 = 50 # kmol/h
nBenzene_stream10 = 0 # kmol/h
nMethane_stream10 = 0 # kmol/h

T_stream11 = 25 +273.15
nHydrogen_stream11 = 220.9 # kmol/h

# Distillation column
T_condenser_in = 121.8 +273.15 # K; from Siva
T_condenser_out = 120.8 +273.15 # K; from VLE data/Siva
T_reboiler = 152.1 +273.15 # K; from VLE data/Siva
refluxratio = 1.47 # from Siva
nBenzene_recycle = refluxratio * nBenzene_stream8

nHydrogen_condenser = 0 # kmol/h
nMethane_condenser = 0 # kmol/h
nBenzene_condenser = nBenzene_recycle + nBenzene_stream8 # kmol/h
nToluene_condenser = 0 # kmol/h
massflow_condenser = MassflowProcessfluid(nHydrogen_condenser, nMethane_condenser, \
                                          nBenzene_condenser, nToluene_condenser)

nHydrogen_reboiler = 0 # kmol/h
nMethane_reboiler = 0 # kmol/h
nBenzene_reboiler = nBenzene_stream7 + nBenzene_recycle # kmol/h
nToluene_reboiler = nToluene_stream7 # kmol/h
massflow_reboiler = MassflowProcessfluid(nHydrogen_reboiler, nMethane_reboiler, \
                                        nBenzene_reboiler, nToluene_reboiler)

# Heater
T_fluegasOut = 300 +273.15 # K

# Utilities
# Cooling water
T_cw_in = 30+273.15 # ingoing temperature of cooling water in K
T_cw_out = 45+273.15 # outgoing temperature of cooling water in K (should be 40-45 K according to Turton)
# Low pressure steam (5 barg)
T_lps = 160+273.15 # K
# Medium pressure steam (10 barg)
T_mps = 184+273.15 # K
# High pressure steam (41 barg)
T_hps = 254+273.15 # K








