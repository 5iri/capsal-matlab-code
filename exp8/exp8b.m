clc; clear; close all;

% Given Cost Coefficients
a = [500 400 200];      % Constant coefficients
b = [5.3 5.5 5.8];      % Linear coefficients
c = [0.004 0.006 0.009];% Quadratic coefficients

% Generator limits (MW)
P_min = [200 150 100];  % Minimum power limits for generators
P_max = [450 350 225];  % Maximum power limits for generators

% Total Power Demand (MW)
Pd = 975;

% Initial guess for lambda
lambda = 6.0;

% Convergence tolerance
tolerance = 1e-5;

% Maximum iterations
max_iter = 1000;

% Iteration counter
iter = 0;

while true
    iter = iter + 1;
    
    % Compute power output for each unit from lambda
    P1 = (lambda - b(1)) / (2 * c(1));
    P2 = (lambda - b(2)) / (2 * c(2));
    P3 = (lambda - b(3)) / (2 * c(3));
    
    % Apply generator limits
    P1 = max(min(P1, P_max(1)), P_min(1)); % P1 must be between P_min(1) and P_max(1)
    P2 = max(min(P2, P_max(2)), P_min(2)); % P2 must be between P_min(2) and P_max(2)
    P3 = max(min(P3, P_max(3)), P_min(3)); % P3 must be between P_min(3) and P_max(3)
    
    % Total generated power
    P_total = P1 + P2 + P3;
    
    % Check for convergence (power balance)
    if abs(P_total - Pd) < tolerance
        break;
    end
    
    % Compute dP/dLambda (sensitivity) based on the current power outputs
    dP_dLambda = (1 / (2 * c(1))) + (1 / (2 * c(2))) + (1 / (2 * c(3)));
    
    % Update lambda properly
    lambda = lambda + (Pd - P_total) / dP_dLambda;
    
    if iter > max_iter
        disp('Lambda iteration did not converge.');
        break;
    end
end

% Display results
disp('-----------------------------------------');
disp('      Economic Dispatch Results');
disp('-----------------------------------------');
fprintf('Lambda (Final): %.4f\n', lambda);
fprintf('Generator 1 Output (P1): %.4f MW\n', P1);
fprintf('Generator 2 Output (P2): %.4f MW\n', P2);
fprintf('Generator 3 Output (P3): %.4f MW\n', P3);
fprintf('Total Generation: %.4f MW\n', P_total);
fprintf('Total Demand: %.4f MW\n', Pd);

% Compute total cost
C1 = a(1) + b(1)*P1 + c(1)*P1^2;
C2 = a(2) + b(2)*P2 + c(2)*P2^2;
C3 = a(3) + b(3)*P3 + c(3)*P3^2;
Total_Cost = C1 + C2 + C3;

fprintf('Total Generation Cost: र%.2f\n', Total_Cost);
fprintf('Number of Iterations: %d\n', iter);
disp('-----------------------------------------');
