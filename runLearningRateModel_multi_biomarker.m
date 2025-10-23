%         MULTI-VARIABLE DYNAMIC SIMULATION VIA LEARNING RATE
% This script simulates the kidney's homeostatic response by monitoring
% multiple biomarkers (Na+, K+, HCO3-) and adjusting its reabsorption
% parameters to restore balance.

clc; clear; close all;

%%  1. DEFINE CONSTANTS, TARGETS, AND MODEL PARAMETERS 

%  Physiological Setpoints (The "Healthy" Target Values) 
healthy_conc.Na = 140; healthy_conc.K = 4.25; healthy_conc.HCO3 = 24;
healthy_conc.Urea = 4.75; healthy_conc.Cl = 101; healthy_conc.Glucose = 5.0;

%  Baseline Kidney Parameters 
baseline_GFR = 105; % Baseline GFR in mL/min
% Order: [Na; K; HCO3; Urea; Cl; Glucose; Water]
base_reab.pct      = [0.65; 0.5707; 0.85; 0.25; 0.60; 1; 0.66];
base_reab.desc     = [0; 0; 0; 0.075; 0; 0; 0.15];
base_reab.asc      = [0.25; 0.1756; 0; 0.075; 0.25; 0; 0];
base_reab.dct      = [0.05; 0.02195; 0.085; 0; 0.075; 0; 0];
base_reab.cort     = [0.02; 0.0439; 0.045; 0.0125; 0.035; 0; 0.075];
base_reab.med      = [0.02; 0.03951; 0.005; 0.1125; 0.015; 0; 0.04];

%  Model Control Parameters 
learning_rate_Na = 0.01;   % Learning rate for sodium and water
learning_rate_K = 0.05;    % K+ regulation is more sensitive, so we give it a slightly higher rate
learning_rate_HCO3 = 0.01; % Bicarbonate adjustments are subtle
TOTAL_BODY_WATER_L = 42;

%%  2. DEFINE THE INITIAL "UNHEALTHY" STATE 
% Let's simulate a state of MILD RENAL INSUFFICIENCY, where the body has
% trouble excreting waste products, leading to high Potassium and Urea.
num_days = 20; % Run for longer to see the full correction
initial_unhealthy_conc = healthy_conc;
initial_unhealthy_conc.Na = 138;     % Sodium might be slightly low
initial_unhealthy_conc.K = 5.5;      % High Potassium (Hyperkalemia)
initial_unhealthy_conc.Urea = 15.0;  % High Urea (Azotemia)
initial_unhealthy_conc.HCO3 = 20;    % Low Bicarbonate (Metabolic Acidosis)


%%  3. PRE-CALCULATE THE HEALTHY STATE FOR COMPARISON 
[healthy_streams, ~] = nephronModel_engine(healthy_conc, baseline_GFR / 60, base_reab);
healthy_NaCl_delivery = healthy_streams(4, 2) + healthy_streams(4, 6);
healthy_daily_loss = healthy_streams(7, 2:8) * 24; % mol/day


%%  4. RUN THE MULTI-DAY LEARNING SIMULATION 

% History arrays
history.plasma_Na = zeros(1, num_days + 1);
history.plasma_K = zeros(1, num_days + 1);
history.plasma_HCO3 = zeros(1, num_days + 1);
history.reab_Na = zeros(1, num_days);
history.reab_H2O = zeros(1, num_days);
history.reab_K = zeros(1, num_days);

% Initialize the simulation
current_plasma_conc = initial_unhealthy_conc;
adjusted_reab = base_reab;
history.plasma_Na(1) = current_plasma_conc.Na;
history.plasma_K(1) = current_plasma_conc.K;
history.plasma_HCO3(1) = current_plasma_conc.HCO3;

fprintf(' STARTING MULTI-VARIABLE LEARNING SIMULATION \n');

for day = 1:num_days
    fprintf('\n Day %d \n', day);
    fprintf('Starting Plasma [Na+]: %.1f, [K+]: %.2f, [HCO3-]: %.1f\n', ...
        current_plasma_conc.Na, current_plasma_conc.K, current_plasma_conc.HCO3);

    %  A. The Multi-Variable Control System (The Learning Logic) 
    
    % == SODIUM AND WATER CONTROL ==
    if current_plasma_conc.Na > healthy_conc.Na + 0.5
        adjusted_reab.dct(1)  = adjusted_reab.dct(1) - learning_rate_Na; % Decrease Na reabsorption
        adjusted_reab.cort(7) = adjusted_reab.cort(7) + learning_rate_Na; % Increase H2O reabsorption
    elseif current_plasma_conc.Na < healthy_conc.Na - 0.5
        adjusted_reab.dct(1)  = adjusted_reab.dct(1) + learning_rate_Na; % Increase Na reabsorption
        adjusted_reab.cort(7) = adjusted_reab.cort(7) - learning_rate_Na; % Decrease H2O reabsorption
    end
    
    % == POTASSIUM CONTROL ==
    if current_plasma_conc.K > healthy_conc.K + 0.1
        % To EXCRETE more K+, we increase SECRETION. We simulate this by
        % aggressively DECREASING the reabsorption rate into negative values.
        adjusted_reab.cort(2) = adjusted_reab.cort(2) - learning_rate_K;
        adjusted_reab.med(2) = adjusted_reab.med(2) - learning_rate_K;
    elseif current_plasma_conc.K < healthy_conc.K - 0.1
        % To CONSERVE more K+, we INCREASE reabsorption.
        adjusted_reab.cort(2) = adjusted_reab.cort(2) + learning_rate_K;
        adjusted_reab.med(2) = adjusted_reab.med(2) + learning_rate_K;
    end

    % == BICARBONATE (pH) CONTROL ==
    if current_plasma_conc.HCO3 < healthy_conc.HCO3 - 0.5 % Acidosis
        % To correct acidosis, body must RETAIN all bicarbonate.
        % We do this by INCREASING the already high PCT reabsorption rate.
        adjusted_reab.pct(3) = adjusted_reab.pct(3) + learning_rate_HCO3;
    elseif current_plasma_conc.HCO3 > healthy_conc.HCO3 + 0.5 % Alkalosis
        % To correct alkalosis, body must EXCRETE some bicarbonate.
        % We do this by slightly DECREASING the PCT reabsorption rate.
        adjusted_reab.pct(3) = adjusted_reab.pct(3) - learning_rate_HCO3;
    end
    
    % Boundary checks to keep reabsorption rates realistic
    fields = fieldnames(adjusted_reab);
    for i = 1:length(fields)
        % For K+, we allow negative reabsorption (secretion), so we give it a different lower bound.
        if strcmp(fields{i}, 'cort') || strcmp(fields{i}, 'med')
            adjusted_reab.(fields{i})(2) = max(-5.0, min(0.99, adjusted_reab.(fields{i})(2)));
        end
        adjusted_reab.(fields{i}) = max(0, min(0.99, adjusted_reab.(fields{i})));
    end

    %  B. Tubuloglomerular Feedback (TGF) Simulation 
    % (This section remains the same, providing GFR stability)
    current_gfr = baseline_GFR;
    for iter = 1:5
        [daily_streams, ~] = nephronModel_engine(current_plasma_conc, current_gfr / 60, adjusted_reab);
        current_NaCl_delivery = daily_streams(4, 2) + daily_streams(4, 6);
        tgf_factor = healthy_NaCl_delivery / current_NaCl_delivery;
        current_gfr = baseline_GFR * tgf_factor;
        current_gfr = max(90, min(120, current_gfr));
    end
    final_daily_gfr = current_gfr;
    
    %  C. UPDATE PLASMA FOR THE NEXT DAY 
    [final_streams, ~] = nephronModel_engine(current_plasma_conc, final_daily_gfr / 60, adjusted_reab);
    unhealthy_daily_loss = final_streams(7, 2:8) * 24;
    net_daily_change_moles = healthy_daily_loss - unhealthy_daily_loss;
    total_body_moles = (TOTAL_BODY_WATER_L * [current_plasma_conc.Na, current_plasma_conc.K, ...
        current_plasma_conc.HCO3, current_plasma_conc.Urea, current_plasma_conc.Cl, current_plasma_conc.Glucose] / 1000);
    total_body_water_moles = TOTAL_BODY_WATER_L * 1000 / 18;
    new_total_body_moles = total_body_moles + net_daily_change_moles(1:6);
    new_total_body_water_moles = total_body_water_moles + net_daily_change_moles(7);
    new_total_body_water_L = new_total_body_water_moles * 18 / 1000;
    new_plasma_conc_mmol_L = (new_total_body_moles ./ new_total_body_water_L) * 1000;
    current_plasma_conc.Na = new_plasma_conc_mmol_L(1);
    current_plasma_conc.K = new_plasma_conc_mmol_L(2);
    current_plasma_conc.HCO3 = new_plasma_conc_mmol_L(3);

    %  D. Store Daily History 
    history.plasma_Na(day+1) = current_plasma_conc.Na;
    history.plasma_K(day+1) = current_plasma_conc.K;
    history.plasma_HCO3(day+1) = current_plasma_conc.HCO3;
    history.reab_Na(day) = adjusted_reab.dct(1);
    history.reab_H2O(day) = adjusted_reab.cort(7);
    history.reab_K(day) = adjusted_reab.cort(2); % Track Cortical CD K+ rate
end


%%  5. PLOT THE MULTI-DAY RESULTS 
days_vector = 1:num_days;
days_vector_plasma = 0:num_days;

figure('Name', 'Multi-Variable Homeostasis Simulation');

%  Plot 1: NA Plasma Concentrations Over Time 
subplot(2,2,1);
plot(days_vector_plasma, history.plasma_Na, '-s', 'LineWidth', 2, 'DisplayName', 'Plasma [Na+]');
hold on;
% Plot healthy setpoints for reference
plot(days_vector_plasma, ones(1, num_days+1) * healthy_conc.Na, 'k--');
title('NA Plasma Concentrations vs. Time');
ylabel('Concentration (mmol/L)');
grid on; legend('show', 'Location', 'best');

%  Plot 2: K Plasma Concentrations Over Time 
subplot(2,2,2);
hold on;
plot(days_vector_plasma, history.plasma_K, '-^', 'LineWidth', 2, 'DisplayName', 'Plasma [K+]');
% Plot healthy setpoints for reference
plot(days_vector_plasma, ones(1, num_days+1) * healthy_conc.K, 'k--');
title('K Plasma Concentrations vs. Time');
ylabel('Concentration (mmol/L)');
grid on; legend('show', 'Location', 'best');

%  Plot 3: HCO3 Plasma Concentrations Over Time 
subplot(2,2,3);
hold on;
plot(days_vector_plasma, history.plasma_HCO3, '-d', 'LineWidth', 2, 'DisplayName', 'Plasma [HCO3-]');
% Plot healthy setpoints for reference
plot(days_vector_plasma, ones(1, num_days+1) * healthy_conc.HCO3, 'k--');
title('HCO3 Plasma Concentrations vs. Time');
ylabel('Concentration (mmol/L)');
grid on; legend('show', 'Location', 'best');

%  Plot 4: Learned Reabsorption Rates 
subplot(2,2,4);
plot(days_vector, history.reab_Na * 100, '-o', 'LineWidth', 2, 'DisplayName', 'DCT Na+ Reabsorption (%)');
hold on;
plot(days_vector, history.reab_H2O * 100, '-^', 'LineWidth', 2, 'DisplayName', 'CD H2O Reabsorption (%)');
plot(days_vector, history.reab_K * 100, '-d', 'LineWidth', 2, 'DisplayName', 'CD K+ Reabsorption/Secretion (%)');
title('Learned Reabsorption Rates Over Time');
ylabel('Reabsorption / Secretion (%)');
xlabel('Day');
grid on; legend('show', 'Location', 'best');