% ---------------------------------------------------------------
% Practical 2 – Step-by-step formulation of Y_BUS  (Direct Method)
% ---------------------------------------------------------------
% Line data: [From  To   R(p.u.)  X(p.u.)]
ld = [ ...
    1  2   0.01   0.15 ;  % line 1-2
    1  3   0.02   0.25 ;  % line 1-3
    1  4   0.03   0.35 ;  % line 1-4
    2  3   0.03   0.35 ;  % line 2-3
    3  4   0.01   0.15 ;  % line 3-4
    4  5   0.04   0.50 ]; % line 4-5

nbus  = max(ld(:,1:2),[],'all');      % total number of buses
Ybus  = complex(zeros(nbus));         % initialise to 0 + j0

for k = 1:size(ld,1)
    f = ld(k,1);                      % from-bus
    t = ld(k,2);                      % to-bus
    R = ld(k,3);  X = ld(k,4);        % line impedance
    z = R + 1i*X;                     % Z = R + jX
    y = 1/z;                          % admittance

    % Off-diagonal elements (mutual)
    Ybus(f,t) = Ybus(f,t) - y;
    Ybus(t,f) = Ybus(f,t);            % symmetric

    % Diagonal elements (self)
    Ybus(f,f) = Ybus(f,f) + y;
    Ybus(t,t) = Ybus(t,t) + y;
end

disp('Y_BUS (pu):');
disp(Ybus);
