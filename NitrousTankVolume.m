clear; clc;
burntime = 3.96;
T = 90 % [F]
P = 800 % [psi]


m_required = 2.5788 * burntime; % [kg]
P_tank = convpres(P, 'psi', 'pa');
T_tank = convtemp(T, 'f','k'); % [K] N2O temperature
ullage_frac = 0.20;
residual_frac = 0.05;

fluid = 'NitrousOxide';
PropsSI = @(varargin) double(py.CoolProp.CoolProp.PropsSI(varargin{:}));

% Critical point check
T_crit = PropsSI('Tcrit', fluid);
P_crit = PropsSI('Pcrit', fluid);
if T_tank >= T_crit
    error('T = %.2f K is above N2O critical temperature (%.2f K): no liquid phase.', ...
        T_tank, T_crit);
end

% Saturation pressure at tank temperature
P_sat = PropsSI('P', 'T', T_tank, 'Q', 0, fluid)

% Liquid density
if P_tank > P_sat
    rho_liq = PropsSI('D', 'P', P_tank, 'T', T_tank, fluid);
    state   = 'subcooled liquid';

else
    rho_liq = PropsSI('D', 'T', T_tank, 'Q', 0, fluid);
    state   = sprintf(['P < Psat -> using saturated liquid at T ' ...
        '(tank self-pressurizes to %.2f bar)'], P_sat/1e5);
    warning('Specified P (%.2f psi) is below Psat(T) = %.2f psi.', ...
        convpres(P_tank, 'pa','psi'), convpres(P_sat,'pa','psi'));
end



% Sizing
m_residual = m_required * residual_frac;
m_loaded   = m_required + m_residual;
V_liquid   = m_loaded / rho_liq;
V_tank     = V_liquid / (1 - ullage_frac); 
V_ullage   = V_tank - V_liquid;

m3_to_L   = 1e3;
m3_to_in3 = 61023.7441;

fprintf('===== N2O Tank Sizing =====\n');
fprintf('Required usable mass : %8.2f kg\n', m_required);
fprintf('Psat(T)              : %8.2f psia  \n', convpres(P_sat, 'pa','psi'));
fprintf('Critical point       : %8.2f F, %.2f psi\n',...
    convtemp(T_crit,'k','f'), convpres(P_crit,'pa','psi'));
fprintf('State                : %s\n', state);
fprintf('---------------------------\n');
fprintf('Liquid density       : %8.2f kg/m^3\n', rho_liq);
fprintf('Total mass           : %8.3f kg\n', m_loaded);
fprintf('Liquid volume        : %8.3f L\n', V_liquid*m3_to_L);
fprintf('Ullage volume (%2.0f%%)  : %8.3f L\n', ullage_frac*100, V_ullage*m3_to_L);
fprintf('Tank Volume          : %8.3f L \n', V_tank*m3_to_L);

ID = convlength(5.75, "in", "m"); % [in] tank inner diameter
A= pi/4 * ID^2; % [m^2] internal cross-section

L_tank = convlength(V_tank/A, 'm', 'in'); % [m] tank length
h_liquid = convlength(V_liquid / A, 'm', 'in'); % [m] liquid height
h_ullage = convlength(V_ullage / A, 'm', 'in'); % [m] gas length

ID_in = 5.75; % [in] tank inner diameter
ID    = convlength(ID_in, 'in', 'm');  % [m]
A     = pi/4 * ID^2;

L_tank   = convlength(V_tank   / A, 'm', 'in');   % [in] internal length
h_liquid = convlength(V_liquid / A, 'm', 'in');   % [in] liquid height
h_ullage = convlength(V_ullage / A, 'm', 'in');   % [in] gas height

fprintf('---------------------------\n');
fprintf('Total length : %8.2f in \n', L_tank);
fprintf('Liquid length : %8.2f in (%2.0f%%)\n', h_liquid, 100*h_liquid/L_tank);
fprintf('Ullage length : %8.2f in (%2.0f%%)\n', h_ullage, 100*h_ullage/L_tank);
