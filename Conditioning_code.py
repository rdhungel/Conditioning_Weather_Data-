# # This algorithm is for adjusting ambient weather Data to reference conditions

# Name:         Conditioning_Code.py
# Purpose:      Conditioning various desert meterological variables to match agricultural environment
# Based on:		Procedure for "Adjusting Ambient Weather Data to Reference Conditions" by Richard G. Allen, University  of Idaho (White Paper)
# Author:       Ramesh Dhungel - (2010-2011)
# Modfied by:   Bibha Dhungana (2016-10)
# Created       2016-12-16
# Python:       2.7
# Comments:		The code can be be effective probably using functions and redundancy (needs to verfied). 
				# This algorithm is benchmarked in Idaho (spreadsheet) and Colorado from ASCE Manual.

# ---------------------------------------------------------------------------------------------------------------------------
import os

import matplotlib.pyplot as plt
import xlrd
import pandas as pd
import numpy as np
from math import pi
import pylab

# Input from csv file
df = pd.read_csv('C:/Conditioning_Code/Salt desert_aug_metric_2011_2013_test.csv')

####################################################################################################################
#Constants
# Surface emissivity
Emissivity = 0.97

# Specific heat capacity of moist air
Cp = 1013

# Acceleration due to gravity (m/s2)
g = 9.81

# Stefan Bolzman Constant (MJK-4m-2d-1)
Sigma = 0.000000057

# Elevation (m)
DEM = 1580.4

# Landuse type
Landuse = 82

# von Karmon Constant
k = 0.41

# Zero plance displacement (m)
d = 0.05

# Roughness length of momentum (m)
zom = 0.02

# Roughness length of heat (m)
zoh = 0.002

## Height of measurement (m)
z_m = 2

# Height of wind measurement (m)
z_u = 3

# Specific gas constant J / (kg * K)
R = 287

# Blending height 50m
zbl = 50

# Time step
time_step = 1.0

# Lz
Lz = 105

# Latitude in radian
lat_rad = 0.681

# Longitude in degree
long_deg = 114.057

# Albedo
Albedo = 0.25

#####################################################################################################
# Variables based on White Paper (variables names in the code and white paper may differ)
# Ta_amb. = ambient near surface air temperature at zm height (oC)
# ea_amb. = ambient near surface vapor pressure at zm height (kpa)
# Rs	= solar radiation (Watts/m2)
# uz_amb.	= wind speed at zu height (m/s)
# Rl_out	= emitted long-wave radiation under ambient conditions, estimated from Ts_amb (Watts/m2).
# Rl_in	= downwelling long-wave radiation under ambient conditions, estimated from Ts_amb (Watts/m2).
# Ts_amb	= estimated surface temperature under ambient conditions, extrapolated from Ta amb. using estimated surface energy balance components (oC).
# Hamb.	= estimated sensible heat flux under ambient conditions(Watts/m2)
# Rn_amb.	= estimated net radiation under ambient conditions(Watts/m2)
# Gamb.	= estimated soil heat flux under ambient conditions (Watts/m2)
# Tbl	= estimated air temperature at the blending height (oC)
# ebl	= estimated vapor pressure at the blending height (kPa)
# ubl	= estimated wind speed at the blending height (m/s)
# Tref	= estimated air temperature at measurement height under reference conditions (oC)
# eref	= estimated vapor pressure at measurement height under reference conditions (kPa)
# uz ref	= estimated wind speed at measurement height under reference conditions (m/s)
########################################################################################################################

# Temperature is K
df['TKmax'] = df['Ta_amb'] + 273.16

# Solar declination angle
df['decl_rad'] =  0.409 * np.sin ( 2 * pi / 365 * df['DOY'] - 1.39 )

# Inverse relative distance Earth-Sun (correction for eccentricity of Earth's orbit around the sun)
df['dr'] = 1 + 0.033 * np.cos ( 2 * pi * df['DOY'] / 365 )

df['b'] = 2 * pi * ( df['DOY'] - 81 ) / 364

# seasonal correction for solar time [hour].
df ['Sc'] = 0.1645 * np.sin ( 2 * df['b'] ) - 0.1255 * np.cos ( df['b'] ) - 0.025 * np.sin (df['b'])

# Omega s (radians)
df['omega_s'] = np.arccos ( -np.tan (lat_rad) * np.tan (df['decl_rad'] ))

# Omega at 2pm (radians)
df['omega_2pm'] = (pi / 12) * ((( 14 - ( time_step / 2 )) + 0.06667 * ( Lz - long_deg )  + df['Sc'] ) - 12 )

# Solar time angle 1/2 hour before omega, that is, at the beginning of period (radians)
df['omega_1'] = df['omega_2pm'] - pi * time_step / 24

# Solar time angle 1/2 hour after  omega, at end of period (radians)
df['omega_2'] = df['omega_2pm'] + pi * time_step / 24

df['j'] = df['omega_s'] * -1

# Limit of omega 1 (radians)
df['omega_1_lim_max'] =(df[['omega_1', 'j']]).max(axis=1)

df['omega_1_lim'] =(df[['omega_1_lim_max', 'omega_s']]).min(axis=1)

df['k'] = df['omega_s'] * -1

# Limit of omega 2 (radians)
df['omega_2_lim_max'] =(df[['omega_2', 'k']]).max(axis=1)

df['omega_2_lim'] =(df[['omega_2_lim_max', 'omega_s']]).min(axis=1)

# Extraterrestial solar radiation at 2pm (W/m2)
df['Ra'] =  (12 / pi) * 4.92 * df['dr'] * ((df['omega_2'] - df['omega_1']) * np.sin (lat_rad) * np.sin (df['decl_rad']) + np.cos (lat_rad) * np.cos (df['decl_rad']) * (np.sin (df['omega_2_lim']) - np.sin (df['omega_1_lim']))) / 3600 * 1000000

#  Ratio of incoming shortwave radiation to extra terrestial radiation(Rs / Ra)
df['tow_sw'] = df['Rs'] /df['Ra']

# Pressure (kPa)
Pressure = 101.3 * np.power ((( 293 - 0.0065 * DEM ) / 293 ), 5.26 )

# Effective emissivity
df['Eeff'] = 0.75 * np.power (( - 1 * np.log (df['tow_sw'] - 0.001 )) , 0.2 )

# Outgoing longwave radiation (W/m2)
df['Rl_out'] = Emissivity * Sigma * np.power (( df['TKmax'] ) , 4 )

# Incoming longwave radiation(W/m2)
df['Rl_in'] = df['Eeff'] * Sigma * np.power (( df['TKmax'] ) , 4 )

# Estimated net radiation under ambient conditions(W/m2)
df['Rn_amb'] = ( 1 - Albedo ) * df['Rs'] - df['Rl_out'] + df['Rl_in']

#  Intial sensible heat flux
df['Hini'] = df['Rn_amb'] * 0

# Actual Et (W/m2)
df['ETact'] = 0.05 * df['Rn_amb']

# Ambient ground heat flux (W/m2)
df['Gamb'] = ( 0.1 + 0.3 * ( df['Hini'] / df['Rn_amb'] )) * df['Rn_amb']

# Density of air (kg/m)
df['Rho_air'] = ( 1000 * Pressure ) / ( R * 1.01 * ( df['Ta_avg_C'] + 273.15 ))

# Latent heat of vaporization (J / kg)
df['Lambda'] = ( 2.501 - 0.00236 * (df['TKmax']- 273.15 )) * np.power(( 10 ), 6 )

# Actual ET (mm/ sec)
df['ETact_mm_sec'] = df['ETact'] / df['Lambda']

# Ambient sensible heat flux (W/m2)
df['Hamb'] = df['Rn_amb'] - df['Gamb'] - df['ETact']

# Friction velocity netural
df['u*_neu'] = ( k * df['uz_amb'] ) / np.log ( ( z_u - d )/ zom )

# Monin Obukhov Length (m)
df['Monin'] = (df['Rho_air'] * np.power(( - 1 * df['u*_neu'] ),3)) / (k * g * (( df['Hamb'] / ( df['TKmax'] * Cp )) + 0.61 * df['ETact_mm_sec'] ))

# Monin Obukhov parameter
df['one'] =  df['u*_neu']/df['u*_neu']
df['ten'] =  df['one'] * 10

df['x2m_1'] = np.power((1 - ((16 * (2 - 0.05)) / df['Monin'])), 0.25)
df['x2m_2'] =  (df[['Monin', 'one']]).max(axis=1)
df['x2m'] = np.where(df['Monin'] < 0,df['x2m_1'], df['x2m_2'])

# Correction of momentum
df['shim2m_1'] = 2 * np.log (( 1 + df['x2m'] ) / 2 ) + np.log (( 1 + np.power((df['x2m'] ) , 2 )) / 2 ) - ( 2 * np.arctan ( df['x2m'] )) + ( 0.5 * pi )
df['shim2m_2'] = - 5 * ( 2 / df['Monin'] )
df['shim2m'] = np.where(df['Monin'] < 0,df['shim2m_1'], df['shim2m_2'])

# Correction of heat
df['shih2m_1'] = 2 * np.log (( 1 + np.power (( df['x2m'] ) , 2 )) / 2 )
df['shih2m_2'] = - 5 * ( 2 /df['Monin'] )
df['shih2m'] = np.where(df['Monin'] < 0,df['shih2m_1'],df['shih2m_2'] )

# Corrected friction velocity (m/s)
df['u*'] = ( 0.41 * df['uz_amb'] ) / ( np.log ( ( 3 - 0.05 ) / zom ) - df['shim2m'] )

# Corrected aerodynamic resistance (s/m)
df['rah'] = ((np.log (( 2 - 0.05 ) / zoh ) ) - df['shih2m'] ) / ( 0.41 * df['u*'] )

iteration = 1
Not_converged = True

# while Not_converged:
while (iteration < 6):

    df['u*_new'] = df['u*']
    df['rah_new'] = df['rah']

    df['Ts_amb'] = ((df['Hamb'] * df['rah_new']) / (df['Rho_air'] * Cp)) + df['TKmax']
    df['Rl_out'] = Emissivity * Sigma * np.power((df['Ts_amb']), 4)
    df['Rl_in'] = df['Eeff'] * Sigma * np.power((df['Ts_amb']), 4)

    df['Rn_amb'] = (1 - Albedo) * df['Rs'] - df['Rl_out'] + df['Rl_in']

    df['Gamb'] = (0.1 + 0.3 * (df['Hamb'] / df['Rn_amb'])) * df['Rn_amb']
    df['Hamb'] = df['Rn_amb'] - df['Gamb'] - df['ETact']
    df['Monin'] = (df['Rho_air'] * np.power((- 1 * df['u*']), 3)) / (k * g * ((df['Hamb'] / (df['TKmax'] * Cp)) + 0.61 * df['ETact_mm_sec']))

    df['x2m_1'] = np.power((1 - ((16 * (2 - 0.05)) / df['Monin'])), 0.25)
    df['x2m_2'] = (df[['Monin', 'one']]).max(axis=1)
    df['x2m'] = np.where(df['Monin'] < 0, df['x2m_1'], df['x2m_2'])

    # Correction of momentum
    df['shim2m_1'] = 2 * np.log((1 + df['x2m']) / 2) + np.log((1 + np.power((df['x2m']), 2)) / 2) - (2 * np.arctan(df['x2m'])) + (0.5 * pi)
    df['shim2m_2'] = - 5 * (2 / df['Monin'])
    df['shim2m'] = np.where(df['Monin'] < 0, df['shim2m_1'], df['shim2m_2'])

    # Correction of heat
    df['shih2m_1'] = 2 * np.log((1 + np.power((df['x2m']), 2)) / 2)
    df['shih2m_2'] = - 5 * (2 / df['Monin'])
    df['shih2m'] = np.where(df['Monin'] < 0, df['shih2m_1'], df['shih2m_2'])

    # Corrected friction velocity (m/s)
    df['u*'] = (0.41 * df['uz_amb']) / (np.log((3 - 0.05) / zom) - df['shim2m'])

    # Corrected aerodynamic resistance (s/m)
    df['rah'] = ((np.log((2 - 0.05) / zoh)) - df['shih2m']) / (0.41 * df['u*'])

    df['rah_conv'] = np.absolute(df['rah'] - df['rah_new'])

    if  [df['rah_conv'] >= 0.5]:

        Not_converged = True
        iteration += 1
    else:
        Not_converged = False

# Monin from the iterative procedure above
# For x = 50m
df['x50m_1'] = np.power((1 - ((16 * (50 - 0.05)) / df['Monin'])), 0.25)
df['x50m_2'] = (df[['Monin', 'one']]).max(axis=1)
df['x50m'] = np.where(df['Monin'] < 0, df['x50m_1'], df['x50m_2'])

# Correction of heat for x = 50m
df['shih50m_1'] = 2 * np.log((1 + np.power((df['x50m']), 2)) / 2)
df['shih50m_2'] = - 5 * (2 / df['Monin'])
df['shih50m'] = np.where(df['Monin'] < 0, df['shih50m_1'], df['shih50m_2'])

# For x = 10m
df['x10m_1'] = np.power((1 - ((16 * (10 - 0.469)) / df['Monin'])), 0.25)
df['x10m_2'] = (df[['Monin', 'ten']]).max(axis=1)
df['x10m'] = np.where(df['Monin'] < 0, df['x10m_1'], df['x10m_2'])

# Correction of momentum
df['shim10m_1'] = 2 * np.log((1 + df['x10m']) / 2) + np.log((1 + np.power((df['x10m']), 2)) / 2) - (2 * np.arctan(df['x10m'])) + (0.5 * pi)
df['shim10m_2'] = - 5.2 * (( 10 - 0.469 ) / df['x10m'] )
df['shim10m'] = np.where(df['Monin'] < 0, df['shim10m_1'], df['shim10m_2'])

# Temperature at blending height (oC)
df['Tbl'] = df['Ta_amb'] - (df['Hamb'] * ((np.log (( 50 - d ) / ( z_m - d ))) - df['shih50m'] + df['shih2m'] )) / ( 0.41 * df['u*'] * df['Rho_air'] * 1013 )

# Specific humidity ambient
df['qamb'] = 0.622 * (( df['ea_amb'])/ ( Pressure - 0.378 * df['ea_amb'] ))

# Specific humidity at blending height (kg/kg)
df['qbl'] = df['qamb'] - ( df['ETact_mm_sec'] * ((np.log (( zbl - d ) / ( z_m - d ))) - df['shih50m'] + df['shih2m'] )) / ( 0.41 * df['u*'] * df['Rho_air'] )

# vapor pressure at blending height (kPa)
df['ebl'] = Pressure * df['qbl']  / (0.622 + 0.378 * df['qbl'])

# Wind speed at blending height (m/s)
df['ubl'] = df['uz_amb'] + ( df['u*'] * (( np.log (( zbl - d ) / ( 3 - d ))) - df['shim10m'] + df['shim2m'] )) / k

#####################################################################################################################
# Start of new step
df['Rnl'] = ( 170 - (33 * pow  ((df['ea_amb']), 0.5))) * df['tow_sw']

# Net radiation (W/m2)
df['Rn'] = ( 1 - Albedo ) * df['Rs'] - df['Rnl']

# Net radiation (MJ m^-2 h^-1)
df['Rn_1'] = df['Rn'] * 0.0036

# Ground heat flux at reference (W/m2)
df['Gref'] = 0.1 * df['Rn']

# Ground heat flux at reference (MJ m^-2 h^-1)
df['Gref_1'] = 0.1 * df['Rn_1']

# Aerodynamic resistance (s/m)
df['ra_ref'] =  ( np.log (( z_u - (0.67 * 0.12 ))/( 0.12* 0.12 )))* ( np.log (( z_m - (0.67 * 0.12 ))/( 0.012* 0.12 ))) / ( 0.41 * 0.41 * df['uz_amb'] )

# Slope of saturation vapor pressure
df['delta'] = (2503 * (np.exp((17.27 * df['Ta_amb']) / (df['Ta_amb'] + 237.3)))) / pow((df['Ta_amb'] + 237.3), 2)

# Saturation vapor pressure (kPa)
df['es'] = 0.6108 * (np.exp (( 17.27 * df['Ta_amb'] ) / ( df['Ta_amb']+ 237.3 )))

# Latent heat of vaporization
df['Lambda'] = ( 2.501 - ( 0.00236 * df['Ta_amb'] )) * 1000000

# Psychrometric constant (kPa/oC)
df['gamma'] = 1013 * (( Pressure / df['Lambda'] )) / 0.622

# ETo under ambient conditions(W/m2)
df['ETo_amb'] = ((df['delta'] * ( df['Rn'] - df['Gref'] )) + (df['Rho_air'] * 1013 * (( df['es'] - df['ea_amb'] ) / df['ra_ref'] ))) / ( df['delta'] + ( df['gamma'] * ( 1 + ( 50 / df['ra_ref']))))

# Numerator constant that changes with reference type and calculation time step
Cn = 66

# Denominator constant that changes with reference type and calculation time step
Cd = 0.25

# Reference Evapotranspiration for tall surface taken from ASCE Reference manual ETr(mm/hr)
df['ETr_amb'] = ((0.408*df['delta']*( df['Rn_1']-df['Gref_1'] ))+(  df['gamma'] * ( Cn/(df['Ta_avg_C'] + 273))*df['uz_amb']* (df['es'] - df['ea_amb'])))/( df['delta']+ ( df['gamma']*( 1+(Cd*df['uz_amb']))))

# Sensible heat flux at refernce (W/m2)
df['Href'] = df['Rn'] - df['Gref'] - df['ETo_amb']

# Friction velocity at neutral (m/s)
df['u*_neutral'] = ( 0.41 * df['ubl'] ) / ( np.log (( 50 - 0.05 ) / 0.02 ))

# Monin Obukhov length (m)
df['Monin'] = ( df['Rho_air'] * np.power (( - 1 * df['u*_neutral'] ),3)) / (k * g * (( df['Href'] / (( df['Ta_amb'] + 273.16) * Cp )) + (( 0.61 * df['ETo_amb']) / df['Lambda']) ))

# For x = 10m
df['x10m_1'] = np.power((1 - ((16 * (10 - 0.469)) / df['Monin'])), 0.25)
df['x10m_2'] = (df[['Monin', 'ten']]).max(axis=1)
df['x10m'] = np.where(df['Monin'] < 0, df['x10m_1'], df['x10m_2'])

# Correction of momentum
df['shim10m_1'] = 2 * np.log((1 + df['x10m']) / 2) + np.log((1 + np.power((df['x10m']), 2)) / 2) - (2 * np.arctan(df['x10m'])) + (0.5 * pi)
df['shim10m_2'] = - 5.2 * (( 10 - 0.469 ) / df['x10m'] )
df['shim10m'] = np.where(df['Monin'] < 0, df['shim10m_1'], df['shim10m_2'])

# For x = 50m
df['x50m_1'] = np.power((1 - ((16 * (50 - 0.469)) / df['Monin'])), 0.25)
df['x50m_2'] = (df[['Monin', 'ten']]).max(axis=1)
df['x50m'] = np.where(df['Monin'] < 0, df['x50m_1'], df['x50m_2'])

# Correction of heat for x = 50m
df['shih50m_1'] = 2 * np.log((1 + np.power((df['x50m']), 2)) / 2)
df['shih50m_2'] = - 5.2 * ((50 - 0.469)/ df['x50m'])
df['shih50m'] = np.where(df['Monin'] < 0, df['shih50m_1'], df['shih50m_2'])

df['x2m_1'] = np.power((1 - ((16 * (2 - 0.469)) / df['Monin'])), 0.25)
df['x2m_2'] = (df[['Monin', 'ten']]).max(axis=1)
df['x2m'] = np.where(df['Monin'] < 0, df['x2m_1'], df['x2m_2'])

# Correction of heat
df['shih2m_1'] = 2 * np.log((1 + np.power((df['x2m']), 2)) / 2)
df['shih2m_2'] = - 5.2 * (( 2 - 0.469 ) / df['x2m'] )
df['shih2m'] = np.where(df['Monin'] < 0, df['shih2m_1'], df['shih2m_2'])

# Correction of momentum
df['shim2m_1'] = 2 * np.log((1 + df['x2m']) / 2) + np.log((1 + np.power((df['x2m']), 2)) / 2) - (2 * np.arctan(df['x2m'])) + (0.5 * pi)
df['shim2m_2'] = - 5.2 * (( 2 - 0.469 ) / df['x2m'] )
df['shim2m'] = np.where(df['Monin'] < 0, df['shim2m_1'], df['shim2m_2'])

# specific humidity at reference level
df['qa_ref'] = df['qbl'] + (( df['ETo_amb'] / df['Lambda'] ) * (( np.log (( zbl - 0.469 ) / ( z_m - 0.469 ))) - df['shih50m'] + df['shih2m'] )) / ( 0.41 * df['u*_neutral'] * df['Rho_air'] )

# Vapor pressure at reference level
df['e_ref'] = Pressure * ( df['qa_ref'] / ( 0.622 + ( 0.378 * df['qa_ref'] )))

# Air temperature at reference level
df['Ta_ref'] = df['Tbl'] + ( df['Href'] * ((np.log (( 50 - 0.469 ) / ( z_m - 0.469 ))) - df['shih50m'] + df['shih2m'] )) / ( 0.41 * df['u*_neutral'] * df['Rho_air'] * 1013 )

# Wind speed at reference level
df['uz_ref'] = df['ubl'] - ( df['u*_neutral'] * (( np.log (( zbl - 0.469 ) / ( 3 - 0.469 ))) - df['shim10m'] + df['shim2m'] )) / k

print(' Start of iteration from z to zbl')

iteration = 1

Not_converged = True
while (iteration < 8):
# while Not_converged:

    df['e_ref_new'] = df['e_ref']
    df['Ta_ref_new'] = df['Ta_ref']
    df['ra_ref_new'] = df['ra_ref']

    df['Rnl'] = (170 - (33 * pow((df['e_ref_new']), 0.5))) * df['tow_sw']

    # df['Rn'] = (1 - Albedo) * df['Rs_2pm'] - df['Rnl']
    df['Rn'] = (1 - Albedo) * df['Rs'] - df['Rnl']

    # Net radiation (MJm-2h-1)
    df['Rn_1'] = df['Rn'] * 0.0036

    df['Gref'] = 0.1 * df['Rn']

    # Ground heat flux at reference (MJm-2h-1)
    df['Gref_1'] = 0.1 * df['Rn_1']

    df['ra_ref'] = (np.log((z_u - (0.67 * 0.12)) / (0.12 * 0.12))) * (np.log((z_m - (0.67 * 0.12)) / (0.012 * 0.12))) / (0.41 * 0.41 * df['uz_ref'])

    df['delta'] = (2503 * (np.exp((17.27 * df['Ta_ref']) / (df['Ta_ref'] + 237.3)))) / pow((df['Ta_ref'] + 237.3), 2)

    df['es'] = 0.6108 * (np.exp((17.27 * df['Ta_ref']) / (df['Ta_ref'] + 237.3)))

    df['Lambda'] = (2.501 - (0.00236 * df['Ta_ref_new'])) * 1000000

    df['gamma'] = 1013 * ((Pressure / df['Lambda'])) / 0.622

    # ETo_conditioned (Watts/m2)
    df['ETo_conditioned'] = ((df['delta'] * (df['Rn'] - df['Gref'])) + (df['Rho_air'] * 1013 * ((df['es'] - df['e_ref']) / df['ra_ref_new']))) / (
    df['delta'] + (df['gamma'] * (1 + (50 / df['ra_ref_new']))))

    # ETr_Conditioned(mm/hour)
    df['ETr_conditioned'] = ((0.408 * df['delta'] * (df['Rn_1'] - df['Gref_1'])) + (df['gamma'] * (Cn / (df['Ta_avg_C'] + 273)) * df['uz_amb'] * (df['es'] - df['e_ref']))) / (df['delta'] + (df['gamma'] * (1 + (Cd * df['uz_amb']))))


    df['Href'] = df['Rn'] - df['Gref'] - df['ETo_conditioned']
    df['u*_neutral'] = (0.41 * df['ubl']) / (np.log((50 - 0.469) / 0.084))
    df['Monin'] = (df['Rho_air'] * np.power((- 1 * df['u*']), 3)) / (k * g * ((df['Href'] / ((df['Ta_ref'] + 273.16) * Cp)) + ((0.61 * df['ETo_conditioned']) / df['Lambda'])))

    # For x = 50m
    df['x50m_1'] = np.power((1 - ((16 * (50 - 0.469)) / df['Monin'])), 0.25)
    df['x50m_2'] = (df[['Monin', 'ten']]).max(axis=1)
    df['x50m'] = np.where(df['Monin'] < 0, df['x50m_1'], df['x50m_2'])

    # For x = 10m
    df['x10m_1'] = np.power((1 - ((16 * (10 - 0.469)) / df['Monin'])), 0.25)
    df['x10m_2'] = (df[['Monin', 'ten']]).max(axis=1)
    df['x10m'] = np.where(df['Monin'] < 0, df['x10m_1'], df['x10m_2'])

    #For x = 2m
    df['x2m_1'] = np.power((1 - ((16 * (2 - 0.469)) / df['Monin'])), 0.25)
    df['x2m_2'] = (df[['Monin', 'ten']]).max(axis=1)
    df['x2m'] = np.where(df['Monin'] < 0, df['x2m_1'], df['x2m_2'])

    df['shim10m_1'] = 2 * np.log((1 + df['x10m']) / 2) + np.log((1 + np.power((df['x10m']), 2)) / 2) - (2 * np.arctan(df['x10m'])) + (0.5 * pi)
    df['shim10m_2'] = - 5.2 * ((10 - 0.469) / df['x10m'])
    df['shim10m'] = np.where(df['Monin'] < 0, df['shim10m_1'], df['shim10m_2'])

    df['shih50m_1'] = 2 * np.log((1 + np.power((df['x50m']), 2)) / 2)
    df['shih50m_2'] = - 5.2 * ((50 - 0.469) / df['x50m'])
    df['shih50m'] = np.where(df['Monin'] < 0, df['shih50m_1'], df['shih50m_2'])

    df['shih2m_1'] = 2 * np.log((1 + np.power((df['x2m']), 2)) / 2)
    df['shih2m_2'] = - 5.2 * ((2 - 0.469) / df['x2m'])
    df['shih2m'] = np.where(df['Monin'] < 0, df['shih2m_1'], df['shih2m_2'])

    df['shim2m_1'] = 2 * np.log((1 + df['x2m']) / 2) + np.log((1 + np.power((df['x2m']), 2)) / 2) - (2 * np.arctan(df['x2m'])) + (0.5 * pi)
    df['shim2m_2'] = - 5.2 * ((2 - 0.469) / df['x2m'])
    df['shim2m'] = np.where(df['Monin'] < 0, df['shim2m_1'], df['shim2m_2'])

    df['qa_ref'] = df['qbl'] + ((df['ETo_conditioned'] / df['Lambda']) * ((np.log((zbl - 0.469) / (z_m - 0.469))) - df['shih50m'] + df['shih2m'])) / (0.41 * df['u*'] * df['Rho_air'])

    df['e_ref'] = Pressure * (df['qa_ref'] / (0.622 + (0.378 * df['qa_ref'])))
    df['Ta_ref'] = df['Tbl'] + (df['Href'] * ((np.log((50 - 0.469) / (z_m - 0.469))) - df['shih50m'] + df['shih2m'])) / ( 0.41 * df['u*'] * df['Rho_air'] * 1013)
    df['uz_ref'] = df['ubl'] - (df['u*'] * ((np.log((zbl - 0.469) / (3 - 0.469))) - df['shim10m'] + df['shim2m'])) / k

    df['Ta_ref_conv'] = np.absolute(df['Ta_ref'] - df['Ta_ref_new'])
    df['e_ref_conv'] = np.absolute(df['e_ref'] - df['e_ref_new'])
    df['ra_ref_conv'] = np.absolute(df['ra_ref'] - df['ra_ref_new'])

    df['e_ref_new'] = df['e_ref']
    df['Ta_ref_new'] = df['Ta_ref']
    df['ra_ref_new'] = df['ra_ref']

    iteration += 1

    print (iteration)

df['ETo_amb_mm_hr'] =(df['ETo_amb']/df['Lambda'])*3600

df['ETo_conditioned_mm_hr'] = (df['ETo_conditioned']/df['Lambda'])*3600


df[['Date', 'Time', 'Year', 'DOY', 'e_ref', 'Ta_ref', 'uz_ref', 'qa_ref', 'ETact', 'ETo_amb_mm_hr', 'ETo_conditioned_mm_hr','ETr_amb', 'ETr_conditioned']].to_csv('C:/Conditioning_Code/Conditioning_Results.csv', sep=',', encoding='utf-8')

print (df)
###########################################################################################################################################################
