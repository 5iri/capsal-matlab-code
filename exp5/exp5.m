%% ---------- Data ----------
linedata = [
    1 2 0.07 0.15
    1 3 0.06 0.10
    1 4 0.08 0.25
    2 4 0.04 0.10
    3 4 0.04 0.20
];

nBus  = 4;
slack = 1;
pv    = [2];       % <-- PV bus
pq    = [3 4];     % <-- PQ buses
ang   = [pv pq];   % buses with angle unknowns (non-slack)

% Flat start (or keep V,delta from your NR loop)
V     = [1.05; 1.00; 1.00; 1.00];
delta = zeros(nBus,1);

% Build Ybus and split to G,B
[Ybus, ~, ~] = ybus_calc(linedata, nBus);
G = real(Ybus);  B = imag(Ybus);

%% ---------- Compute Pcalc, Qcalc at current (V,delta) ----------
Pcalc = zeros(nBus,1); Qcalc = zeros(nBus,1);
for i = 1:nBus
    for k = 1:nBus
        d = delta(i) - delta(k);
        Pcalc(i) = Pcalc(i) + V(i)*V(k)*( G(i,k)*cos(d) + B(i,k)*sin(d) );
        Qcalc(i) = Qcalc(i) + V(i)*V(k)*( G(i,k)*sin(d) - B(i,k)*cos(d) );
    end
end

%% ---------- Build Jacobian blocks for PV@2, PQ@3-4 ----------
nang = numel(ang);   % = 3  (δ2, δ3, δ4)
npq  = numel(pq);    % = 2  (V3, V4)

H = zeros(nang, nang);   % ∂P/∂δ   rows: PV+PQ, cols: PV+PQ
N = zeros(nang, npq );   % ∂P/∂V   rows: PV+PQ, cols: PQ
J = zeros(npq,  nang);   % ∂Q/∂δ   rows: PQ,     cols: PV+PQ
L = zeros(npq,  npq );   % ∂Q/∂V   rows: PQ,     cols: PQ

% ----- H & N (rows over ang = [2 3 4]) -----
for ra = 1:nang
    i  = ang(ra);  Vi = V(i);

    % Diagonal terms (use *calculated* Pi,Qi)
    H(ra,ra) = -Qcalc(i) - B(i,i)*Vi^2;

    % N-diagonal only if i is PQ (PV does not have V unknown)
    idx_i_pq = find(pq==i, 1);
    if ~isempty(idx_i_pq)
        N(ra,idx_i_pq) =  Pcalc(i) + G(i,i)*Vi^2;
    end

    % Off-diagonal vs angles (PV+PQ)
    for ca = 1:nang
        if ca==ra, continue; end
        m  = ang(ca); Vm = V(m);
        d  = delta(i) - delta(m);
        Gim = G(i,m); Bim = B(i,m);
        H(ra,ca) = Vi*Vm*( Gim*sin(d) - Bim*cos(d) );
    end

    % Columns over PQ magnitudes (V3,V4)
    for cb = 1:npq
        m  = pq(cb); Vm = V(m);
        d  = delta(i) - delta(m);
        Gim = G(i,m); Bim = B(i,m);
        N(ra,cb) = Vi*Vm*( Gim*cos(d) + Bim*sin(d) );
    end
end

% ----- J & L (rows only for PQ = [3 4]) -----
for rr = 1:npq
    i  = pq(rr);  Vi = V(i);

    % Diagonals
    J(rr, rr + (find(ang==i,1)-rr)) = 0;  % (just keep structure clear)
    J(rr,:) = 0;  % we'll fill full row below
    L(rr,rr) = Qcalc(i) - B(i,i)*Vi^2;

    % Columns over angles (δ2,δ3,δ4)
    for ca = 1:nang
        m  = ang(ca); Vm = V(m);
        d  = delta(i) - delta(m);
        Gim = G(i,m); Bim = B(i,m);
        J(rr,ca) = -Vi*Vm*( Gim*cos(d) + Bim*sin(d) );
    end

    % Columns over PQ magnitudes (V3,V4)
    for cb = 1:npq
        if rr==cb, continue; end
        m  = pq(cb); Vm = V(m);
        d  = delta(i) - delta(m);
        Gim = G(i,m); Bim = B(i,m);
        L(rr,cb) = Vi*Vm*( Gim*sin(d) - Bim*cos(d) );
    end
end

% ----- Final 5x5 Jacobian -----
Jcb = [H N; J L];

fprintf('PV@2, PQ@3-4  ->  Jcb size %dx%d\n', size(Jcb));
disp(Jcb);
