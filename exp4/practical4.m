%PF via outer fixed-point + your gauss_seidel(A,b,x,tol,maxit)
clear; clc;

%--- line data: [from to R X] (pu)
linedata = [
    1 2 0.02   0.04
    1 3 0.01   0.03
    2 3 0.0125 0.025
];

%--- build Ybus
[Ybus, ~, ~] = ybus_calc(linedata, []);
nBus = size(Ybus,1);

%--- slack voltage
V1 = 1.05*exp(1j*0);

%--- given LOADS at PQ buses (positive P,Q = consumption)
S2_load = (2.566 + 1j*1.1102);
S3_load = (1.386 + 1j*0.452);

%convert to NET INJECTIONS Sspec = source - load
Sspec = zeros(nBus,1);
Sspec(2) = -S2_load;
Sspec(3) = -S3_load;

%--- partition Ybus
Ypp = Ybus(2:end, 2:end);   % unknowns: buses 2..n
Yps = Ybus(2:end, 1);       % coupling to slack

%--- outer PF loop (current-injection fixed-point), inner: your gauss_seidel
V = ones(nBus,1); V(1) = V1;          % flat start
tol_pf   = 1e-8;
max_outer= 200;
for k = 1:max_outer
    V_old = V;

    %build current injections I_spec(V) at PQ buses
    v_unknown = V(2:end);
    I_spec = zeros(nBus-1,1);
    for ii = 1:nBus-1
        i_bus = ii + 1;
        I_spec(ii) = conj(Sspec(i_bus)) / conj(V(i_bus));
    end

    %RHS: b = I_spec - Yps * Vslack
    b = I_spec - Yps * V1;

    %solve Ypp * v = b with YOUR gauss_seidel
    x0 = v_unknown;          % start from previous iterate
    [v_unknown, ~, ~] = gauss_seidel(Ypp, b, x0, 1e-10, 1000);

    %update voltages
    V(2:end) = v_unknown;

    %convergence on voltage change
    if max(abs(V - V_old)) < tol_pf
        break
    end
end
iter = k;

%--- results
V2 = V(2); V3 = V(3);
fprintf('v2 = %.6f ∠ %.4f°\n', abs(V2), angle(V2)*180/pi);
fprintf('v3 = %.6f ∠ %.4f°\n', abs(V3), angle(V3)*180/pi);
fprintf('outer iters = %d, max|ΔV| = %.3e\n', iter, max(abs(V - V_old)));
