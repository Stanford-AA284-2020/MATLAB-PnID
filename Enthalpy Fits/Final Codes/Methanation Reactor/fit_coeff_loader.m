%The script generates best fits
clear all; close all;
x_K1 = CO_Methanation_coeff();
x_K3 = CO2_Methanation_coeff();
x_K2 = RWGS_coeff();
x = Cp_coeff();
x_H1 = CO_Methanation_enthalpy_coeff();
x_H3 = CO2_Methanation_enthalpy_coeff();
x_H2 = RWGS_enthalpy_coeff();