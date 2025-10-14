% BIOE 252: Multi-Unit Kidney Reabsorption Model
% This model simulates the molar flow rates of key chemical species through
% a healthy nephron, structured according to the provided Moles
% Conservation Equations.

clc;
clear;
close all;

%% Variables Definition (as per user request)
% n_total(alpha) = total molar flow rate of stream alpha [mol/hr]
% n_i(alpha) = molar flow rate of species i in stream alpha [mol/hr]
% i = {Na+, K+, HCO3-, Urea, Cl-, H2O}
% x_i(alpha) = molar fraction of species i in stream alpha
% V = total volume going into system

%% Initialization and Knowns
% Based on "TABLE OF MARKERS" and standard physiological values.

% Initial Plasma Concentrations (mmol/L)
conc_Na = 140;      % mmol/L
conc_K = 4.25;      % mmol/L (average of 3.5-5.0)
conc_HCO3 = 24;     % mmol/L (average of 22-26)
conc_Urea = 4.75;   % mmol/L (average of 2.5-7.0)
conc_Cl = 101;      % mmol/L (average of 96-106)

% Glomerular Filtration Rate (GFR)
GFR_L_per_min = 0.125; % 125 mL/min is a typical value
GFR = GFR_L_per_min * 60; % Convert to L/hr

% We start the model at Stream 3 (glomerular filtrate) since its
% composition is known and is the input to the tubular system.
% The molar flow rate ni,3 is calculated as Concentration * GFR.
% Units: (mmol/L) * (L/hr) * (1 mol / 1000 mmol) = mol/hr

n_Na_3    = conc_Na * GFR / 1000;
n_K_3     = conc_K * GFR / 1000;
n_HCO3_3  = conc_HCO3 * GFR / 1000;
n_Urea_3  = conc_Urea * GFR / 1000;
n_Cl_3    = conc_Cl * GFR / 1000;
% Water flow: GFR assumes ~99% is water volume. Molar mass of H2O ~ 18 g/mol.
% 1 L water is ~1000g. So flow is (1000 g/L * GFR L/hr) / 18 g/mol.
n_H2O_3   = (1000 * GFR) / 18;

% Total molar flow rate for stream 3
n_total_3 = n_Na_3 + n_K_3 + n_HCO3_3 + n_Urea_3 + n_Cl_3 + n_H2O_3;

% Initialize data storage matrix: 14 streams, 7 columns (n_total, n_Na, n_K, ...)
streams = zeros(14, 7);
streams(3,:) = [n_total_3, n_Na_3, n_K_3, n_HCO3_3, n_Urea_3, n_Cl_3, n_H2O_3];


%% Moles Conservation Equations

% --- Bowman's capsule as the system ---
% Equation: n3 = n4
% This is a direct pass-through of the filtrate.
streams(4,:) = streams(3,:);

% --- PCT as the system ---
% Equation: n4 = n5 + n6
% Stream 5 is reabsorbed to blood. We calculate n6 = n4 - n5.
reab_pct = [ ...
    0.65;   % Na+
    0.65;   % K+
    0.85;   % HCO3- (average of 0.8-0.9)
    0.50;   % Urea
    0.60;   % Cl- (average of 0.5-0.7)
    0.66    % H2O (average of 0.65-0.67)
];
% Calculate molar flow rates for reabsorbed stream 5
n_species_5 = streams(4, 2:7)' .* reab_pct;
n_total_5 = sum(n_species_5);
streams(5, :) = [n_total_5, n_species_5'];

% Calculate molar flow rates for outgoing stream 6
streams(6, :) = streams(4, :) - streams(5, :);

% --- Loop of Henle Descending as the system ---
% Equation: n6 = n7 + n8
% Stream 7 is reabsorbed water. Solutes are mostly impermeable.
n_species_7 = zeros(1,6);
n_species_7(6) = streams(6, 7) * 0.15; % 15% of remaining water reabsorbed
n_total_7 = sum(n_species_7);
streams(7,:) = [n_total_7, n_species_7];

% Calculate stream 8
streams(8,:) = streams(6,:) - streams(7,:);


% --- Loop of Henle thin ascending limb as the system ---
% Equation: n8 = n9 + n10
% Stream 9 is reabsorbed solutes (passive). Impermeable to water.
% Note: The table value for this section is ~0 for all solutes.
% We will assume a small passive reabsorption of NaCl for model completeness.
n_species_9 = zeros(1,6);
n_species_9(1) = streams(8, 2) * 0.05; % Assume 5% of Na+ diffuses out
n_species_9(5) = streams(8, 6) * 0.05; % Assume 5% of Cl- diffuses out
n_total_9 = sum(n_species_9);
streams(9,:) = [n_total_9, n_species_9];

% Calculate stream 10
streams(10,:) = streams(8,:) - streams(9,:);


% --- Loop of Henle thick ascending limb as the system ---
% NOTE: Adjusted equation n10 = n_reabsorbed + n11 for physiological accuracy.
% This segment actively transports solutes out and is impermeable to water.
reab_tal = [ ...
    0.25;   % Na+
    0.20;   % K+
    0;      % HCO3-
    0;      % Urea
    0.50;   % Cl- (Involved in NKCC2 cotransport with Na and K)
    0       % H2O
];
% Calculate reabsorbed species (to interstitial space, not a numbered stream)
n_reabsorbed_species_tal = streams(10, 2:7)' .* reab_tal;
n_reabsorbed_total_tal = sum(n_reabsorbed_species_tal);

% Calculate outgoing stream 11
streams(11, 2:7) = streams(10, 2:7) - n_reabsorbed_species_tal';
streams(11, 1) = sum(streams(11, 2:7));


% --- DCT as the system ---
% Equation: n11 = n12 + n13
% Stream 12 is reabsorbed. Reabsorption is variable and hormone-dependent.
reab_dct = [ ...
    0.05;   % Na+ (5-10% reabsorption)
    0;      % K+ (can be reabsorbed or secreted)
    0;      % HCO3-
    0;      % Urea (DCT is impermeable to urea)
    0.05;   % Cl-
    0.05    % H2O (5-10% of filtered water)
];
n_species_12 = streams(11, 2:7)' .* reab_dct;
n_total_12 = sum(n_species_12);
streams(12, :) = [n_total_12, n_species_12'];

% Calculate stream 13
streams(13, :) = streams(11, :) - streams(12, :);

% --- Collecting Duct as a system ---
% NOTE: Combined Cortical & Medullary. Adjusted equation n13 = n_reabsorbed + n14.
% Stream 14 is the final urine. Reabsorption is hormone-dependent.
reab_cd = [ ...
    0.03;   % Na+ (2-5% reabsorption)
    0;      % K+
    0;      % HCO3-
    0.10;   % Urea (variable reabsorption)
    0.02;   % Cl- (1-2% reabsorption)
    0.08    % H2O (3-10% of filtered water)
];
n_reabsorbed_species_cd = streams(13, 2:7)' .* reab_cd;
n_reabsorbed_total_cd = sum(n_reabsorbed_species_cd);

% Calculate final urine stream 14
streams(14, :) = streams(13, :) - [n_reabsorbed_total_cd, n_reabsorbed_species_cd'];


%% Results Display

% Display Molar Flow Rates (ni)
fprintf('--------------------------------------------------------------------------------------------\n');
fprintf('                               Molar Flow Rates (mol/hr)\n');
fprintf('--------------------------------------------------------------------------------------------\n');
fprintf('Stream\t n_total\t n_Na+\t\t n_K+\t\t n_HCO3-\t n_Urea\t\t n_Cl-\t\t n_H2O\n');
fprintf('--------------------------------------------------------------------------------------------\n');
for i = 3:14 % Displaying streams relevant to the tubule system
    if streams(i,1) ~= 0 % Only display calculated streams
     fprintf('%d\t\t %5.4f\t\t %5.4f\t\t %5.4f\t\t %5.4f\t\t %5.4f\t\t %5.4f\t\t %5.2f\n', ...
            i, streams(i,1), streams(i,2), streams(i,3), streams(i,4), streams(i,5), streams(i,6), streams(i,7));
    end
end
fprintf('--------------------------------------------------------------------------------------------\n\n');


% Calculate and Display Mole Fractions (xi)
fprintf('---------------------------------\n');
fprintf('      Molar Fractions (xi)\n');
fprintf('---------------------------------\n');

% For Filtrate (Stream 4)
x_i_4 = streams(4, 2:7) / streams(4, 1);
fprintf('Stream 4 (Start of Tubule):\n');
fprintf('  x_Na+   = %f\n', x_i_4(1));
fprintf('  x_K+    = %f\n', x_i_4(2));
fprintf('  x_HCO3- = %f\n', x_i_4(3));
fprintf('  x_Urea  = %f\n', x_i_4(4));
fprintf('  x_Cl-   = %f\n', x_i_4(5));
fprintf('  x_H2O   = %f\n\n', x_i_4(6));

% For Urine (Stream 14)
x_i_14 = streams(14, 2:7) / streams(14, 1);
fprintf('Stream 14 (Final Urine):\n');
fprintf('  x_Na+   = %f\n', x_i_14(1));
fprintf('  x_K+    = %f\n', x_i_14(2));
fprintf('  x_HCO3- = %f\n', x_i_14(3));
fprintf('  x_Urea  = %f\n', x_i_14(4));
fprintf('  x_Cl-   = %f\n', x_i_14(5));
fprintf('  x_H2O   = %f\n', x_i_14(6));
fprintf('---------------------------------\n');




%% ========================================================================
%  3. PLOTTING AND VISUALIZATION
%  ========================================================================

% Define the streams that represent the main path of the tubular fluid
tubular_streams_idx = [4, 6, 8, 10, 11, 13, 14];
tubular_streams_labels = {'Bowman''s', 'End PCT', 'End Desc. LOH', 'End Thin Asc.', 'End Thick Asc.', 'End DCT', 'Urine'};
data = streams(tubular_streams_idx, :);

% --- PLOT 1: Molar Flow Rates of Key Species ---
figure('Name', 'Molar Flow Rates of Key Species');
plot(tubular_streams_idx, data(:,7), 'b-o', 'LineWidth', 2, 'DisplayName', 'Water (H2O)');
hold on;
plot(tubular_streams_idx, data(:,2), 'r-s', 'LineWidth', 2, 'DisplayName', 'Sodium (Na+)');
plot(tubular_streams_idx, data(:,6), 'g-^', 'LineWidth', 2, 'DisplayName', 'Chloride (Cl-)');
hold off;
title('Molar Flow Rates of Key Species Along the Nephron');
xlabel('Stream Number');
ylabel('Molar Flow Rate (mol/hr)');
legend('show');
grid on;
xticks(tubular_streams_idx);
xticklabels(tubular_streams_labels);
xtickangle(30);

% --- PLOT 2: Molar Flow Rates of Solutes (Log Scale) ---
figure('Name', 'Molar Flow Rates of Solutes (Log Scale)');
semilogy(tubular_streams_idx, data(:,2), '-s', 'LineWidth', 2, 'DisplayName', 'Na+');
hold on;
semilogy(tubular_streams_idx, data(:,6), '-^', 'LineWidth', 2, 'DisplayName', 'Cl-');
semilogy(tubular_streams_idx, data(:,3), '-d', 'LineWidth', 2, 'DisplayName', 'K+');
semilogy(tubular_streams_idx, data(:,4), '-p', 'LineWidth', 2, 'DisplayName', 'HCO3-');
semilogy(tubular_streams_idx, data(:,5), '-h', 'LineWidth', 2, 'DisplayName', 'Urea');
hold off;
title('Molar Flow Rates of Solutes Along the Nephron (Log Scale)');
xlabel('Stream Number');
ylabel('Molar Flow Rate (mol/hr)');
legend('show', 'Location', 'southwest');
grid on;
xticks(tubular_streams_idx);
xticklabels(tubular_streams_labels);
xtickangle(30);

% --- PLOT 3: Concentration Along the Nephron ---
% Calculate volume in Liters: (moles H2O * 18 g/mol) / (1000 g/L)
volume_L = data(:,7) * 18 / 1000;
% Calculate concentration in mmol/L: (moles solute / L) * 1000 mmol/mol
conc_Na_tubule = (data(:,2) ./ volume_L); % This is in mol/L
conc_Urea_tubule = (data(:,5) ./ volume_L); % This is in mol/L

figure('Name', 'Concentration of Solutes');
plot(tubular_streams_idx, conc_Na_tubule, 'r-s', 'LineWidth', 2, 'DisplayName', 'Na+ Concentration');
hold on;
plot(tubular_streams_idx, conc_Urea_tubule, 'm-p', 'LineWidth', 2, 'DisplayName', 'Urea Concentration');
hold off;
title('Concentration of Solutes Along the Nephron');
xlabel('Stream Number');
ylabel('Concentration (mol/L)');
legend('show');
grid on;
xticks(tubular_streams_idx);
xticklabels(tubular_streams_labels);
xtickangle(30);

% --- PLOT 4: Fluid Composition (Stacked Bar Chart) ---
% Select key points for composition analysis
comp_idx = [4, 6, 11, 14];
comp_labels = {'Filtrate', 'After PCT', 'After LOH', 'Final Urine'};
comp_data = streams(comp_idx, 2:7); % Get species data for these streams

% Calculate mole fractions for the selected streams
mole_fractions = comp_data ./ sum(comp_data, 2);

figure('Name', 'Fluid Composition at Key Points');
bar(mole_fractions, 'stacked');
title('Relative Molar Composition of Tubular Fluid');
ylabel('Mole Fraction (xi)');
xlabel('Location in Nephron');
xticks(1:length(comp_labels));
xticklabels(comp_labels);
legend('Na+', 'K+', 'HCO3-', 'Urea', 'Cl-', 'H2O', 'Location', 'eastoutside');
grid on;