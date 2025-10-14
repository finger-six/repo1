function [streams] = nephronModel(scenarioName, conc_Na, conc_K, conc_HCO3, conc_Urea, conc_Cl)
%NEPHRONMODEL Simulates reabsorption in the kidney for given plasma concentrations.
%   [streams] = nephronModel(scenarioName, conc_Na, conc_K, conc_HCO3, conc_Urea, conc_Cl)
%
%   INPUTS:
%   scenarioName - A string to label the plots (e.g., 'Healthy State').
%   conc_Na, etc. - Initial plasma concentrations in mmol/L.
%
%   OUTPUTS:
%   streams - A 14x7 matrix containing the molar flow rates for all
%             species in every stream of the model.

%% ========================================================================
%  1. MODEL INITIALIZATION AND CALCULATIONS
%  ========================================================================
fprintf('\n\n====================================================================\n');
fprintf('===== Running Simulation for: %s =====\n', upper(scenarioName));
fprintf('====================================================================\n');

%% Initialization and Knowns
GFR_L_per_min = 0.125;
GFR = GFR_L_per_min * 60; % L/hr

n_Na_3    = conc_Na * GFR / 1000;
n_K_3     = conc_K * GFR / 1000;
n_HCO3_3  = conc_HCO3 * GFR / 1000;
n_Urea_3  = conc_Urea * GFR / 1000;
n_Cl_3    = conc_Cl * GFR / 1000;
n_H2O_3   = (1000 * GFR) / 18;

n_total_3 = n_Na_3 + n_K_3 + n_HCO3_3 + n_Urea_3 + n_Cl_3 + n_H2O_3;

streams = zeros(14, 7);
streams(3,:) = [n_total_3, n_Na_3, n_K_3, n_HCO3_3, n_Urea_3, n_Cl_3, n_H2O_3];


%% Moles Conservation Equations
% --- Bowman's capsule ---
streams(4,:) = streams(3,:);

% --- PCT ---
reab_pct = [0.65; 0.65; 0.85; 0.50; 0.60; 0.66];
n_species_5 = streams(4, 2:7)' .* reab_pct;
streams(5, :) = [sum(n_species_5), n_species_5'];
streams(6, :) = streams(4, :) - streams(5, :);

% --- Descending LOH ---
n_species_7 = zeros(1,6);
n_species_7(6) = streams(6, 7) * 0.15;
streams(7,:) = [sum(n_species_7), n_species_7];
streams(8,:) = streams(6,:) - streams(7,:);

% --- Thin Ascending LOH ---
n_species_9 = zeros(1,6);
n_species_9(1) = streams(8, 2) * 0.05;
n_species_9(5) = streams(8, 6) * 0.05;
streams(9,:) = [sum(n_species_9), n_species_9];
streams(10,:) = streams(8,:) - streams(9,:);

% --- Thick Ascending LOH ---
reab_tal = [0.25; 0.20; 0; 0; 0.50; 0];
n_reabsorbed_species_tal = streams(10, 2:7)' .* reab_tal;
streams(11, 2:7) = streams(10, 2:7) - n_reabsorbed_species_tal';
streams(11, 1) = sum(streams(11, 2:7));

% --- DCT ---
reab_dct = [0.05; 0; 0; 0; 0.05; 0.05];
n_species_12 = streams(11, 2:7)' .* reab_dct;
streams(12, :) = [sum(n_species_12), n_species_12'];
streams(13, :) = streams(11, :) - streams(12, :);

% --- Collecting Duct ---
reab_cd = [0.03; 0; 0; 0.10; 0.02; 0.08];
n_reabsorbed_species_cd = streams(13, 2:7)' .* reab_cd;
streams(14, :) = streams(13, :) - [sum(n_reabsorbed_species_cd), n_reabsorbed_species_cd'];


%% ========================================================================
%  2. RESULTS DISPLAY
%  ========================================================================

fprintf('--------------------------------------------------------------------------------------------\n');
fprintf('                               Molar Flow Rates (mol/hr)\n');
fprintf('--------------------------------------------------------------------------------------------\n');
fprintf('Stream\t n_total\t n_Na+\t\t n_K+\t\t n_HCO3-\t n_Urea\t\t n_Cl-\t\t n_H2O\n');
fprintf('--------------------------------------------------------------------------------------------\n');
for i = 3:14
    if streams(i,1) ~= 0
     fprintf('%d\t\t %5.4f\t\t %5.4f\t\t %5.4f\t\t %5.4f\t\t %5.4f\t\t %5.4f\t\t %5.2f\n', ...
            i, streams(i,1), streams(i,2), streams(i,3), streams(i,4), streams(i,5), streams(i,6), streams(i,7));
    end
end
fprintf('--------------------------------------------------------------------------------------------\n');


%% ========================================================================
%  3. PLOTTING AND VISUALIZATION
%  ========================================================================

tubular_streams_idx = [4, 6, 8, 10, 11, 13, 14];
tubular_streams_labels = {'Bowman''s', 'End PCT', 'End Desc. LOH', 'End Thin Asc.', 'End Thick Asc.', 'End DCT', 'Urine'};
data = streams(tubular_streams_idx, :);

% --- PLOT 1: Molar Flow Rates ---
figure('Name', [scenarioName, ': Molar Flow Rates']);
plot(tubular_streams_idx, data(:,7), 'b-o', 'LineWidth', 2, 'DisplayName', 'Water (H2O)');
hold on;
plot(tubular_streams_idx, data(:,2), 'r-s', 'LineWidth', 2, 'DisplayName', 'Sodium (Na+)');
plot(tubular_streams_idx, data(:,6), 'g-^', 'LineWidth', 2, 'DisplayName', 'Chloride (Cl-)');
hold off;
title([scenarioName, ': Molar Flow Rates Along the Nephron']);
xlabel('Nephron Segment'); ylabel('Molar Flow Rate (mol/hr)');
legend('show'); grid on; xticks(tubular_streams_idx); xticklabels(tubular_streams_labels); xtickangle(30);

% --- PLOT 2: Solute Flow Rates (Log Scale) ---
figure('Name', [scenarioName, ': Solute Flow Rates (Log)']);
semilogy(tubular_streams_idx, data(:,2), '-s', 'LineWidth', 2, 'DisplayName', 'Na+');
hold on;
semilogy(tubular_streams_idx, data(:,6), '-^', 'LineWidth', 2, 'DisplayName', 'Cl-');
semilogy(tubular_streams_idx, data(:,3), '-d', 'LineWidth', 2, 'DisplayName', 'K+');
semilogy(tubular_streams_idx, data(:,4), '-p', 'LineWidth', 2, 'DisplayName', 'HCO3-');
semilogy(tubular_streams_idx, data(:,5), '-h', 'LineWidth', 2, 'DisplayName', 'Urea');
hold off;
title([scenarioName, ': Solute Flow Rates (Log Scale)']);
xlabel('Nephron Segment'); ylabel('Molar Flow Rate (mol/hr)');
legend('show', 'Location', 'southwest'); grid on; xticks(tubular_streams_idx); xticklabels(tubular_streams_labels); xtickangle(30);

% --- PLOT 3: Concentration of All Solutes ---
volume_L = data(:,7) * 18 / 1000;
conc_Na_tubule    = (data(:,2) ./ volume_L);
conc_K_tubule     = (data(:,3) ./ volume_L);
conc_HCO3_tubule  = (data(:,4) ./ volume_L);
conc_Urea_tubule  = (data(:,5) ./ volume_L);
conc_Cl_tubule    = (data(:,6) ./ volume_L);

figure('Name', [scenarioName, ': Solute Concentrations']);
plot(tubular_streams_idx, conc_Na_tubule, '-s', 'LineWidth', 2, 'DisplayName', 'Na+');
hold on;
plot(tubular_streams_idx, conc_Cl_tubule, '-^', 'LineWidth', 2, 'DisplayName', 'Cl-');
plot(tubular_streams_idx, conc_Urea_tubule, '-p', 'LineWidth', 2, 'DisplayName', 'Urea');
plot(tubular_streams_idx, conc_K_tubule, '-d', 'LineWidth', 2, 'DisplayName', 'K+');
plot(tubular_streams_idx, conc_HCO3_tubule, '-h', 'LineWidth', 2, 'DisplayName', 'HCO3-');
hold off;
title([scenarioName, ': Solute Concentrations Along the Nephron']);
xlabel('Nephron Segment'); ylabel('Concentration (mol/L)');
legend('show', 'Location', 'northwest'); grid on; xticks(tubular_streams_idx); xticklabels(tubular_streams_labels); xtickangle(30);

% --- PLOT 4: Fluid Composition ---
comp_idx = [4, 6, 11, 14];
comp_labels = {'Filtrate', 'After PCT', 'After LOH', 'Final Urine'};
comp_data = streams(comp_idx, 2:7);
mole_fractions = comp_data ./ sum(comp_data, 2);

figure('Name', [scenarioName, ': Fluid Composition']);
bar(mole_fractions, 'stacked');
title([scenarioName, ': Relative Molar Composition of Fluid']);
ylabel('Mole Fraction (xi)'); xlabel('Location in Nephron');
xticks(1:length(comp_labels)); xticklabels(comp_labels);
legend('Na+', 'K+', 'HCO3-', 'Urea', 'Cl-', 'H2O', 'Location', 'eastoutside'); grid on;

end