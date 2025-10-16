function [streams] = nephronModel3(scenarioName, conc_Na, conc_K, conc_HCO3, conc_Urea, conc_Cl, conc_Glucose)
% Simulates the 6 main tubule segments as
%   distinct units, including separated cortical and medullary ducts.

%  1. set up

%% Initialize Plasma and GFR
GFR_L_per_min = 0.125;
GFR = GFR_L_per_min * 60; % L/hr

% Initialize a 7x8 matrix for the 7 key streams and 8 species (total + 7)
streams = zeros(7, 8);

%  Calculate Stream 1: Fluid Entering the PCT 
streams(1, 2) = conc_Na * GFR / 1000;      % n_Na+
streams(1, 3) = conc_K * GFR / 1000;       % n_K+
streams(1, 4) = conc_HCO3 * GFR / 1000;    % n_HCO3-
streams(1, 5) = conc_Urea * GFR / 1000;    % n_Urea
streams(1, 6) = conc_Cl * GFR / 1000;       % n_Cl-
streams(1, 7) = 0; % Glucose not finished
streams(1, 8) = (1000 * GFR) / 18;         % n_H2O
streams(1, 1) = sum(streams(1, 2:8));      % n_total

% Units

%  1. PCT [Input: Stream 1 -> Output: Stream 2] 
reab_pct = [0.65; 0.65; 0.85; 0.50; 0.60; 0.00; 0.66];
remaining_pct = 1 - reab_pct;
streams(2, 2:8) = streams(1, 2:8) .* remaining_pct';
streams(2, 1) = sum(streams(2, 2:8));

%  2. Descending Limb [Input: Stream 2 -> Output: Stream 3] 
reab_desc = [0; 0; 0; 0.15; 0; 0; 0.15];
remaining_desc = 1 - reab_desc;
streams(3, 2:8) = streams(2, 2:8) .* remaining_desc';
streams(3, 1) = sum(streams(3, 2:8));

%  3. Ascending Limb [Input: Stream 3 -> Output: Stream 4] 
reab_asc = [0.25; 0.20; 0; 0; 0.45; 0; 0];
remaining_asc = 1 - reab_asc;
streams(4, 2:8) = streams(3, 2:8) .* remaining_asc';
streams(4, 1) = sum(streams(4, 2:8));

%  4. DCT [Input: Stream 4 -> Output: Stream 5] 
reab_dct = [0.075; 0; 0.085; 0; 0.075; 0; 0.075];
remaining_dct = 1 - reab_dct;
streams(5, 2:8) = streams(4, 2:8) .* remaining_dct';
streams(5, 1) = sum(streams(5, 2:8));

%  5. Cortical Collecting Duct [Input: Stream 5 -> Output: Stream 6] 
reab_cort_cd = [0.035; 0; 0.045; 0.025; 0.035; 0; 0.075];
remaining_cort_cd = 1 - reab_cort_cd;
streams(6, 2:8) = streams(5, 2:8) .* remaining_cort_cd';
streams(6, 1) = sum(streams(6, 2:8));

%  6. Medullary Collecting Duct [Input: Stream 6 -> Output: Stream 7] 
reab_med_cd = [0.03; 0; 0; 0.225; 0.015; 0; 0.04];
remaining_med_cd = 1 - reab_med_cd;
streams(7, 2:8) = streams(6, 2:8) .* remaining_med_cd';
streams(7, 1) = sum(streams(7, 2:8));


%  3. results
stream_labels = {'1 (PCT In)', '2 (Desc In)', '3 (Asc In)', '4 (DCT In)', '5 (Cort. CD In)', '6 (Med. CD In)', '7 (Final Urine)'};
species_labels_n = {'n_Na+', 'n_K+', 'n_HCO3-', 'n_Urea', 'n_Cl-', 'n_Glucose', 'n_H2O'};
species_labels_c = {'C_Na+', 'C_K+', 'C_HCO3-', 'C_Urea', 'C_Cl-', 'C_Glucose'};

%  Table 1 & 2: Molar Flow Rates and Concentrations 
volume_L = streams(:,8) * 18 / 1000;
concentrations = streams(:, 2:7) ./ volume_L;
fprintf('\n\n TABLE 1: Molar Flow Rates (mol/hr) \n');
fprintf('%-16s', 'Stream');
for j = 1:length(species_labels_n), fprintf('\t%s', species_labels_n{j}); end, fprintf('\n');
for i = 1:7
    fprintf('%-16s', stream_labels{i});
    fprintf('\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t\t%5.2f\n', streams(i, 2:8));
end
fprintf('\n\n TABLE 2: Solute Concentrations (mol/L) \n');
fprintf('%-16s', 'Stream');
for j = 1:length(species_labels_c), fprintf('\t%s', species_labels_c{j}); end, fprintf('\n');
for i = 1:7
    fprintf('%-16s', stream_labels{i});
    fprintf('\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%5.4f\n', concentrations(i, :));
end


%  4. plotting
stream_indices_for_plotting = 1:7;

%  PLOT 1: Solute Flow Rates (Log Scale) 
figure('Name', [scenarioName, ': Solute Flow Rates (Log)']);
semilogy(stream_indices_for_plotting, streams(:,2), '-s', 'LineWidth', 2, 'DisplayName', 'Na+'); hold on;
semilogy(stream_indices_for_plotting, streams(:,6), '-^', 'LineWidth', 2, 'DisplayName', 'Cl-');
semilogy(stream_indices_for_plotting, streams(:,5), '-p', 'LineWidth', 2, 'DisplayName', 'Urea');
semilogy(stream_indices_for_plotting, streams(:,3), '-d', 'LineWidth', 2, 'DisplayName', 'K+');
semilogy(stream_indices_for_plotting, streams(:,4), '-h', 'LineWidth', 2, 'DisplayName', 'HCO3-');
semilogy(stream_indices_for_plotting, streams(:,7), '-x', 'LineWidth', 2, 'DisplayName', 'Glucose (zero)');
hold off;
title([scenarioName, ': Solute Flow Rates (Log Scale)']);
xlabel('Stream Number (Input to Segment)'); ylabel('Molar Flow Rate (mol/hr)');
legend('show', 'Location', 'southwest'); grid on; xticks(stream_indices_for_plotting); xticklabels(stream_labels); xtickangle(45);

%  PLOT 2: Concentration of All Solutes 
figure('Name', [scenarioName, ': Solute Concentrations']);
plot(stream_indices_for_plotting, concentrations(:,1), '-s', 'LineWidth', 2, 'DisplayName', 'Na+'); hold on;
plot(stream_indices_for_plotting, concentrations(:,5), '-^', 'LineWidth', 2, 'DisplayName', 'Cl-');
plot(stream_indices_for_plotting, concentrations(:,4), '-p', 'LineWidth', 2, 'DisplayName', 'Urea');
plot(stream_indices_for_plotting, concentrations(:,2), '-d', 'LineWidth', 2, 'DisplayName', 'K+');
plot(stream_indices_for_plotting, concentrations(:,3), '-h', 'LineWidth', 2, 'DisplayName', 'HCO3-');
plot(stream_indices_for_plotting, concentrations(:,6), '-x', 'LineWidth', 2, 'DisplayName', 'Glucose (zero)');
hold off;
title([scenarioName, ': Solute Concentrations Along the Nephron']);
xlabel('Stream Number (Input to Segment)'); ylabel('Concentration (mol/L)');
legend('show', 'Location', 'best'); grid on; xticks(stream_indices_for_plotting); xticklabels(stream_labels); xtickangle(45);

end



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



