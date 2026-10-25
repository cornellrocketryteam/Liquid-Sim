clear; clc;

burntime = 3.96; % [s]
P_tank = convpres(600, 'psi', 'pa'); 
T_tank = convtemp(90, 'f', 'k');
ullage_frac = 0.20;
residual_frac = 0.05; 
ID_in = 5.75; % [in] tank ID

fluid   = 'Ethanol';
PropsSI = @(varargin) double(py.CoolProp.CoolProp.PropsSI(varargin{:}));


m_required = 0.8512 * burntime;

% Phase check: ethanol must be liquid (P above its vapor pressure)
T_crit = PropsSI('Tcrit', fluid);
P_sat  = PropsSI('P', 'T', T_tank, 'Q', 0, fluid);
if T_tank >= T_crit
    error('T = %.2f K is above ethanol critical temperature (%.2f K).', T_tank, T_crit);
end
if P_tank < P_sat
    error('P (%.3f bar) is below ethanol vapor pressure (%.3f bar): fuel would boil.', ...
        P_tank/1e5, P_sat/1e5);
end

% Liquid density (subcooled / compressed liquid)
rho_liq = PropsSI('D', 'P', P_tank, 'T', T_tank, fluid);

% Sizing
m_residual = m_required * residual_frac;
m_loaded   = m_required + m_residual;
V_liquid   = m_loaded / rho_liq;
V_tank     = V_liquid / (1 - ullage_frac);
V_ullage   = V_tank - V_liquid;

% Length (flat bulkheads)
ID = convlength(ID_in, 'in', 'm');            % [m]
A  = pi/4 * ID^2;                             % [m^2]
L_tank   = convlength(V_tank   / A, 'm', 'in'); % [in] internal length
h_liquid = convlength(V_liquid / A, 'm', 'in'); % [in] liquid height
h_ullage = convlength(V_ullage / A, 'm', 'in'); % [in] gas height

m3_to_L   = 1e3;
m3_to_in3 = 61023.7441;

fprintf('===== Ethanol Tank Sizing =====\n');
fprintf('Required combustion mass : %8.3f kg\n', m_required);
fprintf('Tank pressure        : %8.2f bar  (%.1f psia)\n', P_tank/1e5, P_tank/6894.757);
fprintf('Tank temperature     : %8.2f K    (%.1f F)\n', T_tank, convtemp(T_tank,'k','f'));
fprintf('Vapor pressure       : %8.4f bar  (%.2f psia)\n', P_sat/1e5, P_sat/6894.757);
fprintf('-------------------------------\n');
fprintf('Liquid density       : %8.2f kg/m^3\n', rho_liq);
fprintf('Total loaded mass    : %8.3f kg\n', m_loaded);
fprintf('Liquid volume        : %8.3f L\n', V_liquid*m3_to_L);
fprintf('Ullage volume (%2.0f%%)  : %8.3f L\n', ullage_frac*100, V_ullage*m3_to_L);
fprintf('Tank volume          : %8.3f L  (%.1f in^3)\n', V_tank*m3_to_L, V_tank*m3_to_in3);
fprintf('-------------------------------\n');
fprintf('Internal length      : %8.2f in\n', L_tank);
fprintf('Liquid height        : %8.2f in\n', h_liquid);
fprintf('Ullage height        : %8.2f in\n', h_ullage);