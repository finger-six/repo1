

clc;
clear;
close all;

% Normal starting concentrations.
healthy_Na    = 140;
healthy_K     = 4.25;
healthy_HCO3  = 24;
healthy_Urea  = 4.75;
healthy_Cl    = 101;
[healthy_streams] = nephronModel3('Healthy State', healthy_Na, healthy_K, healthy_HCO3, healthy_Urea, healthy_Cl);
pause; 



