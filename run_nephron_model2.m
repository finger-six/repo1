% =========================================================================
%                    RUNNER SCRIPT FOR NEPHRON MODEL
% =========================================================================
% This script defines and runs different physiological scenarios to test
% the kidney's response using the nephronModel function.

clc;
clear;
close all;

%% --- SCENARIO 1: HEALTHY STATE ---
% Normal physiological concentrations.
disp('RUNNING SCENARIO 1: HEALTHY STATE. Press any key to continue after plots appear...');
healthy_Na    = 140;
healthy_K     = 4.25;
healthy_HCO3  = 24;
healthy_Urea  = 4.75;
healthy_Cl    = 101;
[healthy_streams] = nephronModel('Healthy State', healthy_Na, healthy_K, healthy_HCO3, healthy_Urea, healthy_Cl);
pause; % Waits for you to press a key before running the next scenario


%% --- SCENARIO 2: HIGH SODIUM (HYPERNATREMIA / DEHYDRATION) ---
% Simulates a state of dehydration where sodium concentration is high.
disp('RUNNING SCENARIO 2: HIGH SODIUM. Press any key to continue after plots appear...');
high_Na_Na    = 155;  % High Sodium
high_Na_K     = 4.5;
high_Na_HCO3  = 22;
high_Na_Urea  = 6.0;  % Urea may also be slightly elevated
high_Na_Cl    = 110;
[highNa_streams] = nephronModel('High Sodium', high_Na_Na, high_Na_K, high_Na_HCO3, high_Na_Urea, high_Na_Cl);
pause;


%% --- SCENARIO 3: HIGH UREA (RENAL INSUFFICIENCY) ---
% Simulates a state where the kidneys are not effectively clearing waste,
% leading to a buildup of urea in the blood.
disp('RUNNING SCENARIO 3: HIGH UREA...');
high_Urea_Na    = 140;
high_Urea_K     = 5.0; % K+ might be slightly elevated as well
high_Urea_HCO3  = 24;
high_Urea_Urea  = 25; % Significantly high Urea
high_Urea_Cl    = 101;
[highUrea_streams] = nephronModel('High Urea', high_Urea_Na, high_Urea_K, high_Urea_HCO3, high_Urea_Urea, high_Urea_Cl);


%% --- SCENARIO 4: COMPARATIVE ANALYSIS ---
% This plot compares a key output metric across all simulated scenarios.
disp('PLOTTING COMPARATIVE ANALYSIS...');

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