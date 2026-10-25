% Finds feedforward constants by finding the valve angle provides
% exactly the N2 needed to hold a nitrous/ethanol tank at P_set during the burn
% at several COPV pressures, then least-squares fits theta = theta_0 + b_ff/P (P = COPV pressure [psig]).

% Inputs:
%   P_set_psig - tank setpoint during the burn [psig]
%   mdot_liq - total propellant mass flow out of the tank [kg/s]
%   rho_liq - propellant liquid density [kg/m^3]
%   Pv - vapor pressure [Pa abs]
%   valve_angle - valve angle data [deg]
%   valve_Cv - valve Cv data at valve_angle
%   csv_file - N2 critical flow table
%   P_atm_psi - local atmospheric pressure [psia]
%   P_fit_psig - COPV pressures to fit over [psig]
%   T_fit_K - COPV gas temperatures at those pressures [K]
%   T_ull_K - ullage gas temperature during the burn [K]
%   name - (optional) leg name for the printed report, e.g. 'Ethanol'
% Outputs:
%   theta_0 - feedforward intercept [deg]
%   b_ff - feedforward slope on 1/P [deg*psi]
%   theta_req - required angle at each P_fit_psig [deg]

function [theta_0, b_ff, theta_req] = fitFeedForward(P_set_psig, mdot_liq, rho_liq, Pv, valve_angle, valve_Cv, csv_file, ...
    P_atm_psi, P_fit_psig, T_fit_K, T_ull_K, name)
if nargin < 12, name = ''; end

% Conversions & properties
R = 8.314462618/0.0280134; % N2 gas constant [J/kg-K]
g2a = @(pg) convpres(pg + P_atm_psi,'psi','Pa'); % psig -> Pa absolute
in2m2 = convlength(1,'in','m')^2; % 1 in^2 -> m^2
N2 = @(out,in1,v1,in2,v2) py.CoolProp.CoolProp.PropsSI(out,in1,v1,in2,v2,'Nitrogen');
% C* table
tab = readmatrix(csv_file, 'NumHeaderLines', 0);
Pg = tab(1, 2:end);
Tg = tab(2:end, 1);
Cg = fillmissing(tab(2:end, 2:end), 'linear', 2, 'EndValues', 'nearest');
Cstar = @(T,P) interp2(Pg, Tg, Cg, min(max(P,Pg(1)),Pg(end)), min(max(T,Tg(1)),Tg(end)));

% 1) N2 mass flow and pressure for P_set
P_set = g2a(P_set_psig); % [Pa abs]
P_N2 = max(P_set - Pv, convpres(1,'psi','Pa')); % ullage N2 partial pressure, min 1 psi [Pa abs]
rho_N2_ull = N2('D','P',P_N2,'T',T_ull_K); % ullage N2 density [kg/m^3]
mdot_N2_req = rho_N2_ull * mdot_liq/rho_liq; % [kg/s]

% 2) Required Cv and angle at each COPV state
iOpen = find(valve_Cv == 0, 1, 'last');
theta_req = zeros(size(P_fit_psig));
for i = 1:numel(P_fit_psig)
    P = g2a(P_fit_psig(i)); % [Pa abs]
    Cv_req = 38*mdot_N2_req*sqrt(R*T_fit_K(i)) / (Cstar(T_fit_K(i),P)*P*in2m2);
    if Cv_req > max(valve_Cv)
            error('fitFeedforward (%s): needs Cv = %.3f at %g psig, above the valve max Cv %.3f', name, Cv_req, P_fit_psig(i), max(valve_Cv));
    end
    theta_req(i) = interp1(valve_Cv(iOpen:end), valve_angle(iOpen:end), Cv_req, 'pchip');
end

% Least-squares fit of theta = theta_0 + b_ff*(1/P)
c = polyfit(1./P_fit_psig, theta_req, 1);
b_ff = c(1);
theta_0 = c(2);
% Report
fprintf('=== Feedforward fit %s===\n', [name ' ']);
fprintf('Setpoint %.1f psig, N2 needed %.4f kg/s (ullage %.0f K)\n', P_set_psig, mdot_N2_req, T_ull_K);
for i = 1:numel(P_fit_psig)
    fprintf('  COPV %6.0f psig, %5.1f K: required %.2f deg, fit %.2f deg\n', ...
        P_fit_psig(i), T_fit_K(i), theta_req(i), theta_0 + b_ff/P_fit_psig(i));
end
fprintf('theta_0 = %.2f deg, b_ff = %.0f deg*psi\n\n', theta_0, b_ff);
end