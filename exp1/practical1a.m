% Parameters for the short transmission line
R = 0.01; % Resistance in ohms
X = 0.1; % Reactance in ohms
V_s = 230; % Sending end voltage in volts
P_load = 50; % Load power in watts
V_r = 0.9 * V_s; % Receiving end voltage (90% of sending end voltage)
% Length of the transmission line in kilometers
length_of_line = 50; % Assuming l contains the length in km
power_factor = 0.9; % Power factor
P_load = P_load * power_factor; % Adjust load power based on power factor
% Calculate the current
I = P_load / V_r;

% Calculate voltage drop
V_drop = I * (R + 1j * X);

% Calculate receiving end voltage
V_r_calculated = V_s - V_drop;

% Calculate efficiency
efficiency = (P_load / (P_load + abs(V_drop * I))) * 100;

% Calculate voltage regulation
voltage_regulation = ((V_s - abs(V_r_calculated)) / abs(V_r_calculated)) * 100;

% Display results
fprintf('Calculated Receiving End Voltage: %.4f V\n', abs(V_r_calculated));
fprintf('Efficiency: %.4f%%\n', efficiency);
fprintf('Voltage Regulation: %.4f%%\n', voltage_regulation);