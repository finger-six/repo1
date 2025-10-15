

clc;
clear;
close all;

% Normal physiological concentrations.
healthy_Na    = 140;
healthy_K     = 4.25;
healthy_HCO3  = 24;
healthy_Urea  = 4.75;
healthy_Cl    = 101;
[healthy_streams] = nephronModel2('Healthy State', healthy_Na, healthy_K, healthy_HCO3, healthy_Urea, healthy_Cl);
pause; % Waits for you to press a key before running the next scenario


% Simulates a state of dehydration where sodium concentration is high.
high_Na_Na    = 155;  % High Sodium
high_Na_K     = 4.5;
high_Na_HCO3  = 22;
high_Na_Urea  = 6.0;  % Urea may also be slightly elevated
high_Na_Cl    = 110;
[highNa_streams] = nephronModel2('High Sodium', high_Na_Na, high_Na_K, high_Na_HCO3, high_Na_Urea, high_Na_Cl);
pause;



%% --- SCENARIO 4: COMPARATIVE ANALYSIS ---
% This plot compares a key output metric across all simulated scenarios.

% Get final urine concentration of Urea for each scenario
volume_healthy = healthy_streams(14,7) * 18 / 1000;
urea_conc_healthy = (healthy_streams(14,5) / volume_healthy);

volume_highNa = highNa_streams(14,7) * 18 / 1000;
urea_conc_highNa = (highNa_streams(14,5) / volume_highNa);

volume_highUrea = highUrea_streams(14,7) * 18 / 1000;
urea_conc_highUrea = (highUrea_streams(14,5) / volume_highUrea);

% Data for the comparison plot
comparison_data = [urea_conc_healthy, urea_conc_highNa, urea_conc_highUrea];
scenario_labels = {'Healthy', 'High Sodium', 'High Urea'};

figure('Name', 'Comparative Analysis');
bar(comparison_data);
title('Comparison of Final Urine Urea Concentration');
ylabel('Urea Concentration (mol/L)');
xlabel('Scenario');
xticklabels(scenario_labels);
grid on;