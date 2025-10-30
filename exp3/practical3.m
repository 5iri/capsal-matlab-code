% Power System Analysis Lab - Practical 3
% Formulation of Ybus using Singular Transformation Method
clear; clc;

% ---------------- Line Data ----------------
% Format: [From To R X]
linedata = [
    1 2 0.01 0.15
    1 3 0.02 0.25
    1 4 0.03 0.35
    2 3 0.03 0.35
    3 4 0.01 0.15
    4 5 0.04 0.50
];

nBus = 5;               % number of buses
refBus = 1;             % choose reference bus
m = size(linedata,1);   % number of branches

% -------- Step 1: Branch admittances --------
y = 1 ./ (linedata(:,3) + 1i*linedata(:,4));   % admittance
Yprim = diag(y);                                % primitive admittance

% -------- Step 2: Incidence matrix ----------
A = zeros(m, nBus);
for k = 1:m
    f = linedata(k,1);
    t = linedata(k,2);
    A(k,f) =  1;
    A(k,t) = -1;
end

% Reduced incidence (remove reference bus column)
keep = setdiff(1:nBus, refBus);
K = A(:, keep);

% -------- Step 3: Ybus matrices -------------
Ybus_reduced = K' * Yprim * K;
Ybus_full    = A' * Yprim * A;

% -------- Display results ------------------
disp('Primitive branch admittances:');
disp(diag(Yprim));

disp('Reduced Incidence Matrix K:');
disp(K);

disp('Reduced Ybus (excluding ref bus):');
disp(Ybus_reduced);

disp('Full Ybus:');
disp(Ybus_full);

% -------- Validation -----------------------
fprintf('\nCheck symmetry: %d\n', isequal(round(Ybus_full,10), round(Ybus_full.',10)));
fprintf('Row sums (should be ~0 if no shunts):\n');
disp(sum(Ybus_full,2));
