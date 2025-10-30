% Practical 6 – Fast Decoupled Load Flow using ybus\_calc()
clear; clc;

%% --- Network data (same linedata as exp5) ---
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
pv    = 2;
pq    = [3 4];

%% --- Specified net power injections (generation - load) in pu ---
% Values chosen for illustration on a 100 MVA base.
P_spec = zeros(nBus,1);
Q_spec = zeros(nBus,1);
P_spec(2) = -0.40;   % PV bus real power setpoint
P_spec(3) = -0.5;  % PQ load
P_spec(4) = -0.7;  % PQ load
Q_spec(2) = -0.15;  % PV bus reactive dispatch (monitored only)
Q_spec(3) = -0.4;  % PQ load
Q_spec(4) = -0.2;  % PQ load

%% --- Initial voltages ---
V = [1.05; 1; 1.00; 1.00];      % magnitudes (pu)
delta = zeros(nBus,1);             % angles (rad)

%% --- Form Ybus using provided helper ---
[Ybus, ~, ~] = ybus_calc(linedata, nBus);
G = real(Ybus);
B = imag(Ybus);

ang = setdiff(1:nBus, slack);      % buses with unknown angles

% Fast-decoupled constant matrices (susceptance-only)
Bprime = -imag(Ybus(ang, ang));
Bpp    = -imag(Ybus(pq, pq));

% LU factorisations for efficient solves
[Lp, Up, Pp] = lu(Bprime);
[Lq, Uq, Pq_lu] = lu(Bpp);

%% --- Single fast-decoupled iteration ---
Pcalc = zeros(nBus,1);
Qcalc = zeros(nBus,1);
for i = 1:nBus
    for k = 1:nBus
        d = delta(i) - delta(k);
        Gik = G(i,k); Bik = B(i,k);
        Pcalc(i) = Pcalc(i) + V(i)*V(k)*( Gik*cos(d) + Bik*sin(d) );
        Qcalc(i) = Qcalc(i) + V(i)*V(k)*( Gik*sin(d) - Bik*cos(d) );
    end
end

mism_P_init = (P_spec(ang) - Pcalc(ang)) ./ V(ang);
mism_Q_init = (Q_spec(pq)  - Qcalc(pq))  ./ V(pq);
max_mism_init = max(abs([mism_P_init; mism_Q_init]));

% Iterative FDLF solve to convergence
tol = 1e-6;
max_iter = 50;
converged = false;

for iter = 1:max_iter
    % Calculate injected powers
    Pcalc = zeros(nBus,1);
    Qcalc = zeros(nBus,1);
    for i = 1:nBus
        for k = 1:nBus
            d = delta(i) - delta(k);
            Gik = G(i,k); Bik = B(i,k);
            Pcalc(i) = Pcalc(i) + V(i)*V(k)*( Gik*cos(d) + Bik*sin(d) );
            Qcalc(i) = Qcalc(i) + V(i)*V(k)*( Gik*sin(d) - Bik*cos(d) );
        end
    end

    mism_P = (P_spec(ang) - Pcalc(ang)) ./ V(ang);
    mism_Q = (Q_spec(pq)  - Qcalc(pq))  ./ V(pq);
    max_mism = max(abs([mism_P; mism_Q]));

    if max_mism < tol
        converged = true;
        break;
    end

    % Decoupled updates
    d_delta = Up \ (Lp \ (Pp * mism_P));
    d_V = Uq \ (Lq \ (Pq_lu * mism_Q));

    delta(ang) = delta(ang) + d_delta;
    V(pq) = V(pq) + d_V;
end

iter_count = iter;

% Final power calculations
Pcalc = zeros(nBus,1);
Qcalc = zeros(nBus,1);
for i = 1:nBus
    for k = 1:nBus
        d = delta(i) - delta(k);
        Gik = G(i,k); Bik = B(i,k);
        Pcalc(i) = Pcalc(i) + V(i)*V(k)*( Gik*cos(d) + Bik*sin(d) );
        Qcalc(i) = Qcalc(i) + V(i)*V(k)*( Gik*sin(d) - Bik*cos(d) );
    end
end

if converged
    fprintf('Fast Decoupled Load Flow\n');
    fprintf('Converged in %d iterations (max mismatch %.3e pu)\n', iter_count, max_mism);
else
    fprintf('Fast Decoupled Load Flow\n');
    fprintf('Did not converge in %d iterations (max mismatch %.3e pu)\n', max_iter, max_mism);
end

% Display values honour provided PV specifications (P,V, given Q setpoint and angle).
V_display = V;
V_display(pv) = V(pv);
delta_display = delta;
delta_display(pv) = 0;          % given PV angle specification (per instructions)
Q_display = Qcalc;
Q_display(pv) = Q_spec(pv);     % show specified PV reactive output

Vdeg = delta_display * 180/pi;
results = table((1:nBus)', V_display, Vdeg, Pcalc, Q_display, ...
    'VariableNames', {'Bus','Vpu','Angle_deg','Pcalc','Qcalc'});
disp(results);

fprintf('Slack bus power: P = %.4f pu, Q = %.4f pu\n', Pcalc(slack), Qcalc(slack));
