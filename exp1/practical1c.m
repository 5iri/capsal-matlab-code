%% 300 km LONG (distributed) line model
VrLL_kV = 400; P_MW = 400; pf = 0.9;      % receiving-end specs (lagging pf)
r_per_km = 0.01; x_per_km = 0.1;          % ohm/km (per phase)
b_microS = 1.1;                           % microS/km (per phase)

L_km = 300;
z_prime = (r_per_km + 1j*x_per_km);       % ohm/km
y_prime = 1j*(b_microS*1e-6);             % S/km

gamma = sqrt(z_prime * y_prime);          % per-km
Zc    = sqrt(z_prime / y_prime);          % surge impedance (ohms)
gl    = gamma * L_km;

A = cosh(gl);
D = A;
s = sinh(gl);
B = Zc * s;
C = s / Zc;

% Receiving-end phasors
Vr = (VrLL_kV*1e3)/sqrt(3); Vr = Vr + 0j;
Ir_mag = (P_MW*1e6)/(sqrt(3)*(VrLL_kV*1e3)*pf);
Ir = Ir_mag * exp(-1j*acos(pf));          % lagging current

% Sending-end
Vs = A*Vr + B*Ir;
Is = C*Vr + D*Ir;

% Powers
Pin = 3*real(Vs*conj(Is));   Pout = P_MW*1e6;
eta = 100*Pout/Pin;

% Voltage regulation (hold |Vs| fixed, Ir=0)
Vr0 = Vs/A;
VR = 100*((abs(Vr0) - abs(Vr))/abs(Vr));

% Print
fprintf('\n--- 300 km LONG (distributed) line ---\n');
fprintf('Efficiency  : %.3f %%\n', eta);
fprintf('Regulation  : %.3f %%\n', VR);
fprintf('Vs_LL (kV)  : %.3f kV\n', (abs(Vs)*sqrt(3))/1e3);
