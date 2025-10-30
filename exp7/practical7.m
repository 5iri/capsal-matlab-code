% Practical 7 – DC Power Flow validation using practical 5 system data
clear; clc;

%% --- Network data (same as practical 5) ---
% Columns: from  to  R(pu)  X(pu)
linedata = [
    1 2 0.07 0.15
    1 3 0.06 0.10
    1 4 0.08 0.25
    2 4 0.04 0.10
    3 4 0.04 0.20
];

nBus  = 4;
slack = 1;
ang   = setdiff(1:nBus, slack);  % buses with unknown voltage angles

%% --- Active power injections (generation - load) in pu ---
% Match practical 5 scenario (same injections used in exp6).
P_spec = zeros(nBus,1);
P_spec(2) = -0.40;
P_spec(3) = -0.5;
P_spec(4) = -0.7;
% Sum should be zero so slack injection emerges from solution automatically.

%% --- Build susceptance matrix using existing helper ---
[Ybus, ~, ~] = ybus_calc(linedata, nBus);
B = -imag(Ybus);                  % DC power flow B' matrix (per-unit)
Bred = B(ang, ang);               % reduced matrix excluding slack row/column

%% --- Solve DC power flow ---
theta = zeros(nBus,1);
theta(ang) = Bred \ P_spec(ang);  % solve for non-slack voltage angles (rad)

% Compute resulting power injections from solved angles.
Pcalc = B * theta;
Pslack = Pcalc(slack);

%% --- Line flows (Pij = (theta_i - theta_j)/X_ij) ---
numLines = size(linedata,1);
Pij = zeros(numLines,1);
for k = 1:numLines
    i = linedata(k,1);
    j = linedata(k,2);
    X = linedata(k,4);
    Pij(k) = (theta(i) - theta(j)) / X;
end

%% --- Results ---
fprintf('DC Power Flow (validated on practical 5 system)\n');
fprintf('Slack bus injection: P = %.4f pu\n', Pslack);

results_bus = table(...
    (1:nBus)', ...
    theta * 180/pi, ...
    P_spec, ...
    Pcalc, ...
    'VariableNames', {'Bus','Angle_deg','Pspec','Pcalc'});
disp(results_bus);

line_results = table(...
    linedata(:,1), linedata(:,2), Pij,...
    'VariableNames', {'From','To','Pflow_pu'});
disp(line_results);

% Validation: Pcalc should match specified injections (including slack).
fprintf('Maximum |Pspec - Pcalc| = %.3e pu\n', max(abs(Pcalc - P_spec)));
