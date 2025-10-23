function [streams, concentrations] = nephronModel_engine(input_conc, GFR, reab_rates)
%NEPHRONMODEL_ENGINE Performs a single-pass calculation of the nephron
%   for a given GFR and a given set of reabsorption rates.

%% Initialize Streams based on the provided GFR
streams = zeros(7, 8);
% Stream 1: entering PCT
streams(1, 2) = input_conc.Na * GFR / 1000;
streams(1, 3) = input_conc.K * GFR / 1000;
streams(1, 4) = input_conc.HCO3 * GFR / 1000;
streams(1, 5) = input_conc.Urea * GFR / 1000;
streams(1, 6) = input_conc.Cl * GFR / 1000;
streams(1, 7) = input_conc.Glucose * GFR / 1000; % Glucose is now active
streams(1, 8) = (1000 * GFR) / 18; %h20
streams(1, 1) = sum(streams(1, 2:8));

%% Calculate Remaining Fractions
remaining = {1-reab_rates.pct, 1-reab_rates.desc, 1-reab_rates.asc, ...
             1-reab_rates.dct, 1-reab_rates.cort, 1-reab_rates.med};

%% Compute all downstream flows
for i = 1:6
    streams(i+1, 2:8) = streams(i, 2:8) .* remaining{i}';
    streams(i+1, 1) = sum(streams(i+1, 2:8));
end

%% Compute concentrations
volume_L = streams(:,8) * 18 / 1000;
% Add a small epsilon to prevent division by zero if volume is ever zero
concentrations = streams(:, 2:7) ./ (volume_L + 1e-9);

end