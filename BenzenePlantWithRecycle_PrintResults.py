# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 16:47:40 2023

@author: danfr
"""

from BenzenePlant_functions import *
import BenzenePlantWithRecycle_StreamData as sd
import BenzenePlantWithRecycle_HXEnergyBalances as hx
import BenzenePlantWithRecycle_HeaterEnergyBalances as heater
import BenzenePlantWithRecycle_MechanicalEnergyBalances as pump
import BenzenePlantWithRecycle_CompressorEnergyBalance as comp
import pandas as pd
import openpyxl


# change condenser and reboiler temperature
# check if correct liquid heat enthalpy of toluene/benzene is used
# do pinch analysis

columns = ['E1 (Feed preheater)', 'E2 (Reactor effluent cooler)', 'E3 (Tower feed heater)', \
           'E4 (Benzene condenser)', 'E5 (Product cooler)', 'E6 (Benzene reboiler)', 'H1 (Feed heater)', \
            None, 'E1.1 toluene (liquid)', 'E1.2 hydrogen+methane (vapor)']
index = ['T process stream in [K]', 'T process stream out [K]', 'Mass flow process stream [kg/h]', \
         'Heat flow [MJ/h]', 'UTILITY', 'T utility in [K]', 'T utility out [K]', 'Mass flow cooling water [kg/h]', \
         'Mass flow low pressure steam [kg/h]', 'Mass flow medium pressure steam [kg/h]', \
         'Mass flow high pressure steam [kg/h]', 'Mass flow fuel gas in [kg/h]', 'Mass flow oxygen in [kg/h]', \
         'Mass flow air in [kg/h]', 'Mass flow CO2 out [kh/h]', 'Mass flow water out [kg/h]', None, 'Stream numbers']
T_process_in = [hx.T_stream2[0], sd.T_stream5, sd.T_stream6, sd.T_condenser_in, sd.T_condenser_in, \
                sd.T_reboiler, sd.T_stream3, None, hx.T_stream2liquids[0], hx.T_stream2vapors[0]]
T_process_out = [sd.T_stream3, sd.T_stream6, sd.T_stream7, sd.T_condenser_out, sd.T_stream8, \
                 sd.T_reboiler, sd.T_stream4, None, sd.T_stream3, sd.T_stream3]
m_process = [hx.massflow_stream2liquids+hx.massflow_stream2vapors, sd.massflow_stream5, \
             sd.massflow_stream7, sd.massflow_condenser, sd.massflow_stream8, sd.massflow_reboiler, \
             sd.massflow_stream3, None, hx.massflow_stream2liquids, hx.massflow_stream2vapors]
Q = [hx.Q_E101[0]*10**-3, hx.Q_E102[0]*10**-3, hx.Q_E103*10**-3, hx.Q_E104*10**-3, hx.Q_E105*10**-3, \
     hx.Q_E106*10**-3, heater.Q_H101[0]*10**-3, None, hx.Q_E101liquids*10**-3, hx.Q_E101vapors[0]*10**-3]
T_utility_in = [sd.T_hps, sd.T_cw_in, sd.T_lps, sd.T_cw_in, sd.T_cw_in, sd.T_mps, sd.T_stream6, None, None, None]
T_utility_out = [sd.T_hps, sd.T_cw_out, sd.T_lps, sd.T_cw_out, sd.T_cw_out, sd.T_mps, \
                 sd.T_fluegasOut, None, None, None]
m_cw = [0, hx.massflow_cw_E102, 0, hx.massflow_cw_E104, hx.massflow_cw_E105, 0, 0, None, None, None]
m_lps = [0, 0, hx.massflow_lps_E103, 0, 0, 0, 0, None, None, None]
m_mps = [0, 0, 0, 0, 0, hx.massflow_mps_E106, 0, None, None, None]
m_hps = [hx.massflow_hps_E101, 0, 0, 0, 0, 0, 0, None, None, None]
m_fuelgas = [0, 0, 0, 0, 0, 0, heater.m_fuelgas_furnaceIn, None, None, None]
m_O2 = [0, 0, 0, 0, 0, 0, heater.mO2stoichiometric_furnaceIn, None, None, None]
m_air = [0, 0, 0, 0, 0, 0, heater.mAirstoichiometric_furnaceIn, None, None, None]
m_CO2 = [0, 0, 0, 0, 0, 0, heater.mCO2stoichiometric_furnaceOut, None, None, None]
m_H2O = [0, 0, 0, 0, 0, 0, heater.mH2Ostoichiometric_furnaceOut, None, None, None]
streamnumbers = [2, 5, 7, 'Condenser (C)', 8 , 'Reboiler (R)', None, None, 2, 2]
empty = [None]

data_frame_1 = pd.DataFrame([T_process_in, T_process_out, m_process, Q, empty, \
                             T_utility_in, T_utility_out, m_cw, m_lps, m_mps, m_hps, m_fuelgas, m_O2, \
                             m_air, m_CO2, m_H2O, empty, streamnumbers], index=index, columns=columns)


columns = ['P1', 'P2', 'C1']
index = ['T in [K]', 'T out [K]', 'p in [bar]', 'p out [bar]', 'height [m]', 'Mass flow [kg/h]', 'Shaft work [kW]', 'efficiency']
T_in = [None, None, sd.T_stream6]
T_out = [None, None, comp.T_outCompressor]
p_in = [sd.p_beforeToluenePump*10**-5, 3, sd.p_beforeCompressor*10**-5]
p_out = [sd.p_afterToluenePump*10**-5, 3, sd.p_afterCompressor*10**-5]
height = [0, sd.deltaHeight_RefluxPump, 0]
massflow = [pump.massflow_ToluenePump, pump.massflow_RefluxPump, comp.massflow_Compressor]
Work = [pump.W_P101[0], pump.W_P102, comp.W_Compressor]
efficiency = [sd.efficiency_ToluenePump, sd.efficiency_RefluxPump, sd.efficiency_Compressor]

data_frame_2 = pd.DataFrame([T_in, T_out, p_in, p_out, height, massflow, Work, efficiency],
                  index=index, columns=columns)


with pd.ExcelWriter('EnergyBalances_WithRecycle.xlsx') as writer:
    data_frame_1.to_excel(writer, sheet_name='Heat Exchange with Recycle')
    data_frame_2.to_excel(writer, sheet_name='Shaft Work with Recycle')



