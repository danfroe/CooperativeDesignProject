# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 09:28:01 2023

@author: danfr
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve
from BenzenePlant_functions import *
import BenzenePlantWithRecycle_StreamData as sd


#------------------------------------------------------------------------------
# ENERGY BALANCE
nHydrogen = sd.nHydrogen_stream6-sd.nHydrogen_stream13
nMethane = sd.nMethane_stream6-sd.nMethane_stream13
massflow_Compressor = MassflowProcessfluid(nHydrogen, nMethane, 0, 0)
    
Cp_Compressor = Cp_processfluid_vapor(sd.T_stream6, nHydrogen, nMethane, 0, 0)
var = fsolve(ShaftworkOfCompressor, [100,300], args=(Cp_Compressor, sd.T_stream6, sd.p_beforeCompressor, sd.p_afterCompressor))
w_Compressor = var[0]
T_outCompressor = var[1]
W_Compressor = w_Compressor * massflow_Compressor / 3600 /sd.efficiency_Compressor

print('Temperature out of compressor: %i [K]' % T_outCompressor)
print('Work of compressor: %.2f [kW]' % W_Compressor)
print('Massflow through compressor: %.2f [kg/h]' % massflow_Compressor)