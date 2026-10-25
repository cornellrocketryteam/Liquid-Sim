%% Master Sim
% Phase 1 - Pressurization from EREG into nitrous and ethanol tanks
% Phase 2 - burn. big boom

% COPV: adiabatic, real-gas (CoolProp) blowdown feeding both eth and nitrous lines
% EREG BVs: choked flow from the real-gas C* table; sim errors out if an open valve unchokes
% E-REG: two independent discrete forward euler timestep controllers 
% DP along Lines: Estimates for now
% Ullages: N2 control volume: N2 enthalpy in, heat to the liquid out
% Ethanol: liquid lumped temp (CoolProp)
% Nitrous: subcooled liquid lumped temp (CoolProp). Liquid evaporates into the growing ullage (ullage stays saturated) and cools itself
% Injector: Nitrous(Dyer, CdA set at the design point), ethanol (incompressible orifice), Film cooling ethanol shares the ethanol manifold but bypasses the pintle impingement
% Chamber: Pc = mdot*c*/At with c* fixed at the design point
% Nozzle: isentropic Cf at fixed gamma 

clear; clc; close all;

%% Inputs
P_atm_psi = 14.7; % atmospheric pressure [psia]
burnTime = 3.96; % target burn time for sizing both tanks [s]

% COPV
P_copv_0_psig = 4500; % starting COPV pressure [psig]
T_copv_0_F = 100; % starting COPV temp [F]
V_copv_L = 10; % COPV volume [L]

valve_angle = [0 90*0.1 90*0.2 90*0.3 90*0.4 90*0.5 90*0.6 90*0.7 90*0.8 90*0.9 90];
valve_Cv = [0 0 0 0.01 0.03 0.1 0.22 0.4 0.62 0.9 1.23]; % AVCO HP1100 V15 data

% E-REG controllers: theta_servo = Kp(P)*e + theta_FF(P)
eregOn = true; 
f_ctrl = 50; % controller update rate: PT reading --> servo command [Hz]
t_delay = 1/f_ctrl; % command latency [s]
servo_rate = 319.7; % 319.7 max ball valve rotation rate [deg/s]

theta_init = 0; % valve angle at t = 0 [deg]
theta_min = 0; % BV lower limit[deg]
theta_max_eth = 70; % ethanol valve upper limit [deg] 
theta_max_ox = 70; % nitrous valve upper limit [deg] 

s = 1;
kp_0_eth = 0.05*s; % Kp(P) = kp_0 + a_kp/P [deg/psi]
a_kp_eth = 300*s; % [deg]
kp_0_ox = 0.05*s;
a_kp_ox = 300*s;
theta_crack = 27; % preburn feedforward while pressurizing [deg]

% theta_FF(P) = theta_0 + b_ff/P [deg], P = COPV PT [psig]
% Fit points derived from the COPV starting state (see "E-REG feedforward fit points"):
% pressures span P_copv_0 down to the lowest COPV pressure that keeps both valves choked, temps follow the isentropic blowdown
n_fit = 8; % number of COPV pressures to fit the feedforward over

T_ull_fit_K = 270; % ullage gas temp during burn [K]
t_press_max = 30; % stop if the tanks haven't reached P_set by this time [s]

% Feed lines
dP_COPVtoEth_psi = 50; % COPV -> ethanol tank line loss [psig]
dP_COPVtoOx_psi = 50; % COPV -> nitrous tank line loss [psig] 

% TCA
mdot_annulus = 0.6786; % annulus ethanol mass flow [kg/s]
mdot_film = 0.1726; % film cooling ethanol mass flow [kg/s];
mdot_ox_des = 2.5788; % pintle nitrous mass flow [kg/s]
Pc_psig = 350; % chamber pressure  [psig]

P_injIn_eth_psig = 780; % ethanol annulus inlet pressure at design flow [psig] 
P_injIn_ox_psig = 900; % nitrous injector inlet pressure at design flow [psig]

dP_ethFeed_psi = 200; % ethanol tank -> injector manifold loss at design flow (MAV + lines + regen) [psi]
dP_oxFeed_psi = 75; % nitrous tank -> injector loss at design flow (MAV + lines) [psi]

t_mav = 0.5; % MAV time to fully open [s]
t_hold = 0.1; % hold after both tanks reach trigger pressure before burn [s]
t_ereg_offset = 0; % E-REG activation relative to burn start (MAV open command) [s]: < 0 before, > 0 after. Both tanks still must be triggered for activation to occur
              
ethTriggerOffset = 0; % ethanol burn-trigger offset from setpoint [psi]
nitTriggerOffset = 0; % nitrous burn-trigger offset from setpoint [psi]


% Chamber & nozzle
Dt_mm = 45.98; % throat diameter [mm]
eps_noz = 4.02; % nozzle exit/throat area ratio
Lstar = 0.600; % characteristic length [m]
gam_c = 1.1637; % combustion gas (RPA)
eta_Cf = 0.985; % thrust coefficient efficiency (14 deg exit) https://ntrs.nasa.gov/api/citations/19910019867/downloads/19910019867.pdf

% Ethanol tank - 80 eth, 20 water by mass
w_EtOH = 0.80; % ethanol mass fraction
rho_fuel_ref = 843.49; % blend density at 20 C, 1 atm [kg/m^3] https://nvlpubs.nist.gov/nistpubs/Legacy/circ/nbscircular19.pdf
resid_eth = 0.05; % burn mass fraction
T_eth_0_F = 77; % starting liquid and ullage temp [F]
P_eth_0_psig = 900; % starting tank pressure, sub cooled [psig]
ullFrac_eth = 0.20; % ullage fraction of tank volume


% Nitrous tank
resid_ox = 0.05; % burn mass fraction
T_ox_0_F = 77; % starting liquid and ullage temp [F]
P_ox_0_psig = 1000; % starting nitrous tank pressure [psig]
ullFrac_ox = 0.20; % ullage tank volume

D_tank_in = 5.75; % tank inner diameter [in]
h_ull = 20; % ullage-to-liquid heat transfer coeff [W/m^2-K]

% Solver
dt = 1e-3; % time step [s]
t_max = 10; % max sim time [s]
T_min = 100; % stop if COPV gas cools to this [K] (bottom of C* table)
r_stop = 0.999; % before the burn: stop if P_valve_out/P_copv reaches this (tank can't reach P_set)

scriptPath = fileparts(mfilename('fullpath'));
csv_file = fullfile(scriptPath, 'Reference Tables', 'N2CriticalFlowCSV_Pa.csv');

%% Conversions & functions
R = 8.314462618/0.0280134; % N2 gas constant [J/kg-K]
g0 = 9.80665; % [m/s^2]
g2a = @(pg) convpres(pg + P_atm_psi,'psi','Pa'); % psig -> Pa absolute
psig = @(Pa) convpres(Pa,'Pa','psi') - P_atm_psi; % Pa absolute -> psig
psi2Pa = @(p) convpres(p,'psi','Pa'); % pressure difference psi -> Pa
Pa2psi = @(p) convpres(p,'Pa','psi'); % pressure difference Pa -> psi
N2 = @(out,in1,v1,in2,v2) py.CoolProp.CoolProp.PropsSI(out,in1,v1,in2,v2,'Nitrogen');
EtOH = @(out,in1,v1,in2,v2) py.CoolProp.CoolProp.PropsSI(out,in1,v1,in2,v2,'Ethanol');
NOX = @(out,in1,v1,in2,v2) py.CoolProp.CoolProp.PropsSI(out,in1,v1,in2,v2,'NitrousOxide');


H2O = @(out,in1,v1,in2,v2) py.CoolProp.CoolProp.PropsSI(out,in1,v1,in2,v2,'Water');
M_EtOH = py.CoolProp.CoolProp.PropsSI('M','Ethanol'); % [kg/mol]
M_H2O = py.CoolProp.CoolProp.PropsSI('M','Water'); % [kg/mol]
x_EtOH = (w_EtOH/M_EtOH)/(w_EtOH/M_EtOH + (1-w_EtOH)/M_H2O); % ethanol mole fraction in the liquid (0.610)
rhoIdeal = @(P,T) 1/(w_EtOH/EtOH('D','P',P,'T',T) + (1-w_EtOH)/H2O('D','P',P,'T',T)); % ideal mixing density [kg/m^3]
f_mix = rho_fuel_ref/rhoIdeal(101325, 293.15); % mixing correction
fuelRho = @(P,T) f_mix*rhoIdeal(P,T); % mixture density [kg/m^3]

satEth = satTableFuel(275:1:350, x_EtOH, w_EtOH, M_EtOH, M_H2O, EtOH, H2O); % triple point = 273.16 K


T_copv_0 = convtemp(T_copv_0_F,'F','K'); % [K]
T_eth_0 = convtemp(T_eth_0_F,'F','K'); % [K]
T_ox_0 = convtemp(T_ox_0_F,'F','K'); % [K]
P_atm = convpres(P_atm_psi,'psi','Pa'); % [Pa abs]

% Engine design point
mdot_eth_des = mdot_annulus + mdot_film; % total ethanol out of the tank (annulus + film) [kg/s]
mdot_tot_des = mdot_eth_des + mdot_ox_des; % total propellant flow [kg/s]
OF_des = mdot_ox_des/mdot_eth_des; % O/F (incl. film)
OF_imp_des = mdot_ox_des/mdot_annulus; % O/F at the pintle impingement (annulus ethanol only, no film)
Pc_des = g2a(Pc_psig); % chamber pressure [Pa abs]
dP_inj_eth_des = g2a(P_injIn_eth_psig) - Pc_des; % ethanol injector drop at design flow [Pa]
dP_inj_ox_des = g2a(P_injIn_ox_psig) - Pc_des; % nitrous injector drop at design flow [Pa]
dP_feed_eth_des = psi2Pa(dP_ethFeed_psi); % ethanol feed drop at design flow [Pa]
dP_feed_ox_des = psi2Pa(dP_oxFeed_psi); % nitrous feed drop at design flow [Pa]
P_set_eth = g2a(P_injIn_eth_psig) + dP_feed_eth_des; % tank pressure that gives design flow [Pa abs]
P_set_ox = g2a(P_injIn_ox_psig) + dP_feed_ox_des; % [Pa abs]

% Chamber & nozzle: c* from the design point so design flow gives design Pc
At = pi/4*(Dt_mm/1000)^2; % throat area [m^2]
Ae = eps_noz*At; % exit area [m^2]
cstar_des = Pc_des*At/mdot_tot_des; % effective characteristic velocity c* [m/s]
Gam = sqrt(gam_c)*(2/(gam_c+1))^((gam_c+1)/(2*(gam_c-1))); % Vandenkerckhove function
tau_c = Lstar/(cstar_des*Gam^2); % chamber fill time constant, L*c*/(R*Tc) [s]
Me = fzero(@(M) (1/M)*((2/(gam_c+1))*(1+(gam_c-1)/2*M^2))^((gam_c+1)/(2*(gam_c-1))) - eps_noz, [1.01 10]); % exit Mach
pe_pc = (1+(gam_c-1)/2*Me^2)^(-gam_c/(gam_c-1)); % exit/chamber pressure ratio
Cf_mom = sqrt(2*gam_c^2/(gam_c-1)*(2/(gam_c+1))^((gam_c+1)/(gam_c-1))*(1 - pe_pc^((gam_c-1)/gam_c))); % momentum Cf
thrustFun = @(Pc) max(eta_Cf*Cf_mom*Pc*At + (pe_pc*Pc - P_atm)*Ae, 0); % thrust [N]
F_des = thrustFun(Pc_des);

% Saturation tables (propellant vapor in the ullage): T, Pv, rho_vap, h_fg, cp_liq
satOx = satTable(NOX, 200:0.5:309);

% Nitrous pintle injector inlet : CdA from the design point
T_crit_ox = py.CoolProp.CoolProp.PropsSI('Tcrit','NitrousOxide'); % N2O critical temp [K]
if T_ox_0 >= T_crit_ox - 3
    error('T_ox_0_F = %g F is at/above N2O critical temp (%.1f F): no liquid phase', T_ox_0_F, convtemp(T_crit_ox,'K','F')); 
end
Pv_ox_0 = interpSat(satOx, 2, T_ox_0); % nitrous vapor pressure at fill temp [Pa]
convpres(Pv_ox_0, 'pa', 'psi')
if g2a(P_ox_0_psig) <= Pv_ox_0
    error('P_ox_0_psig must be above the nitrous vapor pressure (%.1f psig at %g F) to start subcooled', psig(Pv_ox_0), T_ox_0_F);
end
rho_ox_inj = NOX('D','P',g2a(P_injIn_ox_psig),'T',T_ox_0); % liquid density at injector inlet [kg/m^3]
s1 = NOX('S','P',g2a(P_injIn_ox_psig),'T',T_ox_0);
h1 = NOX('H','P',g2a(P_injIn_ox_psig),'T',T_ox_0);

P2_scan = linspace(0.2, 1, 81)*Pv_ox_0; 
G_scan = zeros(size(P2_scan));

% equation 16: https://h2tools.org/sites/default/files/ICHS_import/paper_40.pdf
% Finds choking pressure
for j = 1:numel(P2_scan) % HEM mass flux vs back pressure, peak = HEM critical/choked point
    G_scan(j) = NOX('D','P',P2_scan(j),'S',s1)*sqrt(2*max(h1 - NOX('H','P',P2_scan(j),'S',s1), 0));
end
[~, iMax] = max(G_scan);
r_hem = P2_scan(iMax)/Pv_ox_0; % HEM critical back pressure / Pv

% m_dot = Cd*A*G
CdA_ox = mdot_ox_des/oxFlux(g2a(P_injIn_ox_psig), Pc_des, T_ox_0, rho_ox_inj, Pv_ox_0, r_hem, NOX); % pintle effective area [m^2]

% Pintle impinging
rho_eth_inj = fuelRho(g2a(P_injIn_eth_psig), T_eth_0); % injector inlet density [kg/m^3]
v_ann_des = sqrt(2*dP_inj_eth_des/rho_eth_inj); % annulus ethanol velocity [m/s]
v_pin_des = mdot_ox_des/(CdA_ox*rho_ox_inj); % pintle nitrous velocity at liquid density [m/s]
TMR_des = (mdot_ox_des*v_pin_des)/(mdot_annulus*v_ann_des); % total momentum ratio, radial (pintle) / axial (annulus)

% Tank sizing from required propellant mass
rho_eth_0 = fuelRho(g2a(P_eth_0_psig), T_eth_0); % fuel density [kg/m^3]
rho_ox_0 = NOX('D','P',g2a(P_ox_0_psig),'T',T_ox_0); % nitrous subcooled liquid density at starting pressure [kg/m^3]
rho_ox_set = NOX('D','P',P_set_ox,'T',T_ox_0); % nitrous density at burn-start pressure [kg/m^3]
m_eth_0 = mdot_eth_des*burnTime*(1 + resid_eth); % loaded fuel mass [kg]
m_ox_0 = mdot_ox_des*burnTime*(1 + resid_ox); % loaded nitrous mass [kg]
A_surf = pi/4*convlength(D_tank_in,'in','m')^2; % liquid surface area [m^2]

% E-REG feedforward fit points, derived from the COPV starting state
s_copv_0 = N2('S','P',g2a(P_copv_0_psig),'T',T_copv_0); % fill-state entropy [J/kg-K]
P_dn_max = max(P_set_eth + psi2Pa(dP_COPVtoEth_psi), P_set_ox + psi2Pa(dP_COPVtoOx_psi));
P_fit_lo = P_dn_max/0.5283; % first guess: ideal-gas r_crit

% Guess the COPV pressure, look up gam for N2 at that pressure, see what
% choked pressure that gamma predicts, iterate.
for it = 1:50 % r_crit depends on gamma, which depends on r_crt
    T_lo = N2('T','P',P_fit_lo,'S',s_copv_0);
    gam_lo = N2('Cpmass','P',P_fit_lo,'T',T_lo)/N2('Cvmass','P',P_fit_lo,'T',T_lo);
    P_new = P_dn_max/(2/(gam_lo+1))^(gam_lo/(gam_lo-1));
    done = abs(P_new - P_fit_lo) < 100; % converged to 100 Pa
    P_fit_lo = P_new;
    if done
        break; 
    end
end
if P_fit_lo >= g2a(P_copv_0_psig)
    error('COPV start pressure %.0f psig is below the %.0f psig needed to keep the E-REG valves choked at burn-start tank pressure', ...
        P_copv_0_psig, psig(P_fit_lo));
end
P_fit_psig = linspace(P_copv_0_psig, psig(P_fit_lo), n_fit); % COPV pressures to fit over [psig]
T_fit_K = zeros(1, n_fit);
for j = 1:n_fit
    T_fit_K(j) = N2('T','P',g2a(P_fit_psig(j)),'S',s_copv_0); % isentropic COPV gas temp [K]
end

% E-REG feedforward constants (fitFeedForward.m)
Pv_eth_0 = interpSat(satEth, 2, T_eth_0);
[theta_0_eth, b_ff_eth, thReq_eth] = fitFeedForward(psig(P_set_eth), mdot_eth_des, rho_eth_0, Pv_eth_0, valve_angle, valve_Cv, ...
    csv_file, P_atm_psi, P_fit_psig, T_fit_K, T_ull_fit_K, 'Ethanol');
[theta_0_ox, b_ff_ox, thReq_ox] = fitFeedForward(psig(P_set_ox), mdot_ox_des, rho_ox_set, Pv_ox_0, valve_angle, valve_Cv, ...
    csv_file, P_atm_psi, P_fit_psig, T_fit_K, T_ull_fit_K, 'Nitrous');

if any(thReq_eth > theta_max_eth)
    warning('Ethanol leg needs %.1f deg, above theta_max_eth', max(thReq_eth)); 
end

if any(thReq_ox > theta_max_ox)
    warning('Nitrous leg needs %.1f deg, above theta_max_ox', max(thReq_ox));
end

V_copv = V_copv_L/1000; % COPV volume [m^3]
in2m2 = (convlength(1,'in','m'))^2; % 1 in^2 -> m^2 (CdA = Cv/38 in^2)
N_ctrl = round(1/(f_ctrl*dt)); % solver steps per controller update
iCv0 = find(valve_Cv == 0, 1, 'last'); % last zero-Cv angle: start of the invertible (strictly increasing) part of the Cv curve

% C* table: row 1 = pressures [Pa], column 1 = temperatures [K]
tab = readmatrix(csv_file, 'NumHeaderLines', 0);
Pg = tab(1, 2:end);
Tg = tab(2:end, 1);
Cg = fillmissing(tab(2:end, 2:end), 'linear', 2, 'EndValues', 'nearest');
Cstar = @(T,P) interp2(Pg, Tg, Cg, min(max(P,Pg(1)),Pg(end)), min(max(T,Tg(1)),Tg(end)));

%% Starting state
% COPV
P_copv = g2a(P_copv_0_psig); % [Pa abs]
T_copv = T_copv_0; % [K]
m_copv = N2('D','P',P_copv,'T',T_copv) * V_copv;
U_copv = m_copv * N2('U','P',P_copv,'T',T_copv);
m_copv_0 = m_copv;

% Tanks: tk(1) = ethanol, tk(2) = nitrous
legName = {'Ethanol','Nitrous'};
rhoFun = {fuelRho, @(P,T) NOX('D','P',P,'T',T)}; % liquid density property handles (P[Pa abs], T[K])
sat = {satEth, satOx};
T_0 = [T_eth_0, T_ox_0];
m_0 = [m_eth_0, m_ox_0];
m_resid = [m_eth_0 - mdot_eth_des*burnTime, m_ox_0 - mdot_ox_des*burnTime]; % unusable residual [kg] -> burnout
rho_0 = [rho_eth_0, rho_ox_0];
ullFrac = [ullFrac_eth, ullFrac_ox];
P_0 = [g2a(P_eth_0_psig), g2a(P_ox_0_psig)]; % starting tank pressure [Pa abs]: N2 partial = P_0 - Pv
P_set = [P_set_eth, P_set_ox];
dP_trig = psi2Pa([ethTriggerOffset, nitTriggerOffset]); % burn-trigger offsets [Pa]
P_trig = P_set + dP_trig; % tank pressure that arms the burn [Pa abs]
P_pre = P_set + max(dP_trig, 0); % E-REG target before the burn: raised to the trigger if the offset is positive [Pa abs]
theta_max = [theta_max_eth, theta_max_ox];
kp_0 = [kp_0_eth, kp_0_ox]; a_kp = [a_kp_eth, a_kp_ox];
theta_0 = [theta_0_eth, theta_0_ox]; b_ff = [b_ff_eth, b_ff_ox];
dP_line = psi2Pa([dP_COPVtoEth_psi, dP_COPVtoOx_psi]); % N2 line loss [Pa]

for i = 1:2
    tk(i).m_liq = m_0(i); % liquid mass [kg]
    tk(i).T_liq = T_0(i); % [K]
    tk(i).rho = rho_0(i); % liquid density [kg/m^3]    
    V_liq = m_0(i)/rho_0(i); % liquid volume [m^3]
    tk(i).V = V_liq/(1 - ullFrac(i)); % tank volume [m^3]
    tk(i).V_ull = tk(i).V - V_liq; % ullage volume [m^3]
    tk(i).V_ull_0 = tk(i).V_ull;
    tk(i).Pv = interpSat(sat{i}, 2, T_0(i)); % vapor pressure [Pa]
    tk(i).P = P_0(i); % total tank pressure [Pa abs]
    tk(i).P_N2 = tk(i).P - tk(i).Pv; % N2 partial pressure [Pa abs]
    tk(i).T_ull = T_0(i); % [K]
    tk(i).m_N2 = N2('D','P',tk(i).P_N2,'T',T_0(i)) * tk(i).V_ull; % ullage N2 mass [kg]
    tk(i).U_N2 = tk(i).m_N2 * N2('U','P',tk(i).P_N2,'T',T_0(i)); % ullage N2 energy [J]
    tk(i).m_N2_0 = tk(i).m_N2;
    tk(i).m_vap = interpSat(sat{i}, 3, T_0(i)) * tk(i).V_ull; % propellant vapor in ullage [kg]
    tk(i).mdot = 0; % propellant flow out [kg/s]

    % E-REG controller state
    tk(i).theta = theta_init; % actual valve angle [deg]
    tk(i).theta_target = theta_init; % angle the servo is moving to [deg]
    tk(i).cmd_queue = zeros(0,2); % For informing when to send next command: [t + t_delay, theta_cmd]
    tk(i).e = NaN; tk(i).Kp = NaN; tk(i).theta_FF = NaN; tk(i).theta_cmd = theta_init;
    tk(i).t_reach = NaN; % time the tank first reaches P_set
end

% Chamber
Pc = P_atm; % [Pa abs]

%% Solver
n = round(t_max/dt);
[t_log, P_copv_log, T_copv_log, m_copv_log, Pc_log, F_log, OF_log, OFimp_log, TMR_log] = deal(nan(n,1));

[mdot_N2_log, r_log, P_vout_log, theta_log, thetaCmd_log, thetaFF_log, e_log, Kp_log, Cv_log, ...
    P_tank_log, Pv_log, T_ull_log, T_liq_log, m_N2_log, m_vap_log, m_liq_log, V_ull_log, ...
    mdot_liq_log, P_injIn_log, dP_inj_log] = deal(nan(n,2));

stopReason = 'reached t_max';
t_burn = NaN; % time the burn starts
burning = false;
t_ereg = NaN; % time the E-REG switches to burn mode (burn start + t_ereg_offset) [s]


for k = 1:n
    t = (k-1)*dt;

    % tank is within the trigger pressure (setpoint + offset) for first time
    for i = 1:2
        if isnan(tk(i).t_reach) && tk(i).P >= P_trig(i) - psi2Pa(1)
            tk(i).t_reach = t;
        end
    end

    % burn activation
    if ~burning && all(~isnan([tk.t_reach])) && t >= max([tk.t_reach]) + t_hold
        burning = true;
        t_burn = t;
    end

    % E-REG burn-mode activation time: scheduled once both tanks arm (burn start is known from then on)
    if isnan(t_ereg) && all(~isnan([tk.t_reach]))
        t_arm = max([tk.t_reach]); % arming time [s]
        t_ereg = t_arm + t_hold + t_ereg_offset;
        if t_ereg < t_arm
            warning('t_ereg_offset = %g s is earlier than arming (-t_hold = %g s): E-REG burn mode starts at arming', t_ereg_offset, -t_hold);
            t_ereg = t_arm;
        end
    end
    eregBurn = ~isnan(t_ereg) && t >= t_ereg - dt/2; % E-REG in burn mode

    % E-REG controllers: run every 1/f_ctrl on the PT readings
    if mod(k-1, N_ctrl) == 0
        P_copv_meas = max(psig(P_copv), 1); 
        for i = 1:2
            if burning
                P_target = P_set(i); 
            else
                P_target = P_pre(i); 
            end 

            tk(i).e = psig(P_target) - psig(tk(i).P); % pressure error [psi]
            tk(i).Kp = kp_0(i) + a_kp(i)/P_copv_meas; % [deg/psi]

            armed = all(~isnan([tk.t_reach])); % both tanks at trigger

            % Feedforward switches to burn mode at t_ereg (burn start + t_ereg_offset) and ramps over t_mav:
            % N2 demand follows liquid outflow, so the burn feedforward's Cv follows the MAV opening profile, shifted by the offset
            if eregBurn
                f_ff = min((t - t_ereg)/max(t_mav, eps), 1); % feedforward ramp fraction (MAV opening profile, shifted)
                th_burn = min(max(theta_0(i) + b_ff(i)/P_copv_meas, 0), 90); % full-burn feedforward angle [deg]
                Cv_ff = f_ff*interp1(valve_angle, valve_Cv, th_burn, 'pchip'); % feedforward Cv for the current outflow
                tk(i).theta_FF = interp1(valve_Cv(iCv0:end), valve_angle(iCv0:end), Cv_ff, 'pchip'); % back to angle [deg]
            else
                tk(i).theta_FF = theta_crack; % before t_ereg: pressure-hold mode
            end

            tk(i).theta_cmd = min(max(tk(i).Kp*tk(i).e + tk(i).theta_FF, theta_min), theta_max(i));
            if ~eregBurn && tk(i).e <= 0 % before t_ereg: at/above target, hold shut
                tk(i).theta_cmd = theta_min;
            end
            if armed && ~eregOn
                tk(i).theta_FF = NaN; 
                tk(i).theta_cmd = theta_min; 
            end 

            tk(i).cmd_queue(end+1,:) = [t + t_delay, tk(i).theta_cmd];
        end
    end


    % Servo actuation --> limited by servo max angular velocity
    for i = 1:2
        while ~isempty(tk(i).cmd_queue) && tk(i).cmd_queue(1,1) <= t + dt/2
            tk(i).theta_target = tk(i).cmd_queue(1,2);
            tk(i).cmd_queue(1,:) = [];
        end
        tk(i).theta = tk(i).theta + max(min(tk(i).theta_target - tk(i).theta, servo_rate*dt), -servo_rate*dt);
        tk(i).Cv = max(interp1(valve_angle, valve_Cv, tk(i).theta, 'pchip'), 0);
    end

    
    % Upstream (COPV) properties
    cp = N2('Cpmass','P',P_copv,'T',T_copv);
    cv = N2('Cvmass','P',P_copv,'T',T_copv);
    gamma = cp/cv;
    r_crit = (2/(gamma+1))^(gamma/(gamma-1)); % critical pressure ratio
    g_choked = Cstar(T_copv,P_copv) * P_copv / sqrt(R*T_copv); % choked mass flux

    % N2 mass flow through each valve at its current angle
    for i = 1:2
        tk(i).P_valve_out = tk(i).P + dP_line(i); % downstream pressure pressure = tank + line loss [Pa abs]
        tk(i).r = tk(i).P_valve_out/P_copv; % pressure ratio across ereg mav (absolute)
        if tk(i).Cv > 0 && tk(i).r > r_crit 
            error('%s N2 valve UNCHOKED at t = %.3f s: COPV %.1f psig, valve outlet %.1f psig, tank %.1f psig (P_valve_out/P_copv = %.3f > r_crit = %.3f)', ...
                legName{i}, t, psig(P_copv), psig(tk(i).P_valve_out), psig(tk(i).P), tk(i).r, r_crit);
        end
        CdA = tk(i).Cv/38 * in2m2; % effective area [m^2]
        tk(i).mdot_N2 = g_choked * CdA;
    end


    % Engine
    if burning % MAV opening fraction
        if t_mav > 0
            f_mav = min((t - t_burn)/t_mav, 1); 
        else
            f_mav = 1;
        end 

        % Ethanol (annulus + film): P_eth - Pc = dP_feed + dP_inj, both scaling with mdot^2
        if tk(1).P > Pc && f_mav > 0
            tk(1).mdot = mdot_eth_des * sqrt((tk(1).P - Pc)/(dP_feed_eth_des + dP_inj_eth_des/f_mav^2)); %Incompressible orifice equation
        else
            tk(1).mdot = 0;
        end
        tk(1).P_injIn = tk(1).P - dP_feed_eth_des*(tk(1).mdot/mdot_eth_des)^2; % injector inlet pressure [Pa abs]


        % Nitrous (pintle): feed loss from last step's flow, injector flow from the Dyer flux
        tk(2).P_injIn = tk(2).P - dP_feed_ox_des*(tk(2).mdot/mdot_ox_des)^2*(rho_ox_set/tk(2).rho);
        if tk(2).P_injIn > Pc
            tk(2).mdot = f_mav*CdA_ox*oxFlux(tk(2).P_injIn, Pc, tk(2).T_liq, tk(2).rho, tk(2).Pv, r_hem, NOX);
        else
            tk(2).mdot = 0;
        end
        
        % Pintle impingement: radial nitrous out of the pintle vs axial ethanol sheet out of the annulus (film excluded)
        mdot_ann = tk(1).mdot*mdot_annulus/mdot_eth_des; % annulus share of the ethanol flow [kg/s]
        v_ann = sqrt(2*max(tk(1).P_injIn - Pc, 0)/tk(1).rho); % annulus ethanol velocity [m/s]
        v_pin = tk(2).mdot/(max(f_mav, eps)*CdA_ox*tk(2).rho); % pintle nitrous velocity via mass continuity at liquid density [m/s]
        TMR = (tk(2).mdot*v_pin)/(mdot_ann*v_ann); % total momentum ratio

        % Chamber: quasi-steady Pc = mdot*c*/At, approached with the L* fill lag (backward Euler)
        Pc_now = Pc; % Pc the flows above were solved against
        Pc_qs = max((tk(1).mdot + tk(2).mdot)*cstar_des/At, P_atm);
        Pc = Pc + (Pc_qs - Pc)*dt/(tau_c + dt);
        F = thrustFun(Pc_now);

    else
        tk(1).P_injIn = NaN; tk(2).P_injIn = NaN;
        Pc_now = Pc;
        F = 0;
        TMR = NaN;
    end

    % Enthalpy leaving the COPV (and entering the ullages): h = u + P/rho
    % Valve and line losses are adiabatic throttling, so h is unchanged
    h_copv = U_copv/m_copv + P_copv*V_copv/m_copv;

    % Log
    t_log(k) = t;
    P_copv_log(k) = P_copv; T_copv_log(k) = T_copv; m_copv_log(k) = m_copv;
    Pc_log(k) = Pc_now; F_log(k) = F; OF_log(k) = tk(2).mdot/tk(1).mdot;
    OFimp_log(k) = tk(2).mdot/(tk(1).mdot*mdot_annulus/mdot_eth_des); TMR_log(k) = TMR;
    for i = 1:2
        mdot_N2_log(k,i) = tk(i).mdot_N2; r_log(k,i) = tk(i).r; P_vout_log(k,i) = tk(i).P_valve_out;
        theta_log(k,i) = tk(i).theta; thetaCmd_log(k,i) = tk(i).theta_cmd; thetaFF_log(k,i) = tk(i).theta_FF;
        e_log(k,i) = tk(i).e; Kp_log(k,i) = tk(i).Kp; Cv_log(k,i) = tk(i).Cv;
        P_tank_log(k,i) = tk(i).P; Pv_log(k,i) = tk(i).Pv; T_ull_log(k,i) = tk(i).T_ull; T_liq_log(k,i) = tk(i).T_liq;
        m_N2_log(k,i) = tk(i).m_N2; m_vap_log(k,i) = tk(i).m_vap; m_liq_log(k,i) = tk(i).m_liq; V_ull_log(k,i) = tk(i).V_ull;
        mdot_liq_log(k,i) = tk(i).mdot; P_injIn_log(k,i) = tk(i).P_injIn; dP_inj_log(k,i) = tk(i).P_injIn - Pc_now;
    end

    % Stop conditions
    if T_copv < T_min
        stopReason = 'COPV gas below C* table';
        break; 
    end
    if burning && tk(1).m_liq <= m_resid(1)
        stopReason = 'burnout (ethanol down to residual)';
        break; 
    end
    if burning && tk(2).m_liq <= m_resid(2)
        stopReason = 'burnout (nitrous down to residual)';
        break; 
    end
    if burning && any([tk.P] <= Pc_now)
        stopReason = 'tank pressure fell to chamber pressure'; 
        break; 
    end
    if ~burning && any(isnan([tk.t_reach]) & [tk.r] >= r_stop)
        stopReason = 'COPV cannot reach burn-start pressure'; 
        break; 
    end
    if ~burning && t >= t_press_max, stopReason = 'tanks did not reach burn-start pressure by t_press_max';
        break; 
    end

    % COPV: mass + adiabatic energy balance
    mdot_N2_tot = tk(1).mdot_N2 + tk(2).mdot_N2;
    m_copv = m_copv - mdot_N2_tot*dt;
    U_copv = U_copv - mdot_N2_tot*h_copv*dt;

    for i = 1:2
        % Heat from ullage gas into liquid surface
        Q = h_ull * A_surf * (tk(i).T_ull - tk(i).T_liq); % [W]

        % Liquid out through injector --> ullage volume growth this step
        m_liq = tk(i).m_liq - tk(i).mdot*dt;
        V_ull_new = tk(i).V - m_liq/tk(i).rho;
        dV = V_ull_new - tk(i).V_ull; % [m^3]

        % Ullage N2: mass + energy balance (N2 enthalpy in, heat out, boundary work on liquid)

        tk(i).U_N2 = tk(i).U_N2 - Q*dt - tk(i).P_N2*dV + tk(i).mdot_N2*h_copv*dt;
        tk(i).m_N2 = tk(i).m_N2 + tk(i).mdot_N2*dt;

        % Propellant vapor kept at saturation at the liquid temp (instant evaporation) latent heat comes out of the liquid.
        % Linearized-implicit in T_liq (explicit is unstable once the ullage is large)
        T_l = tk(i).T_liq;
        a = interpSat(sat{i}, 3, T_l)*V_ull_new - tk(i).m_vap; % vapor needed at the current temp [kg]
        b = (interpSat(sat{i}, 3, T_l + 0.25) - interpSat(sat{i}, 3, T_l - 0.25))/0.5*V_ull_new; % d(vapor needed)/dT [kg/K]
        h_fg = interpSat(sat{i}, 4, T_l); % [J/kg]
        mc = m_liq * interpSat(sat{i}, 5, T_l); % liquid heat capacity [J/K]
        dT = (Q*dt - h_fg*a)/(mc + h_fg*b); % liquid temp change [K]
        m_evap = a + b*dT; % liquid evaporated this step [kg]
        tk(i).m_vap = tk(i).m_vap + m_evap;
        tk(i).m_liq = m_liq - m_evap;
        tk(i).T_liq = T_l + dT;
        tk(i).V_ull = tk(i).V - tk(i).m_liq/tk(i).rho;

        % New thermodynamic states: N2 from (density, internal energy), vapor from the liquid temp
        tk(i).P_N2 = N2('P','D',tk(i).m_N2/tk(i).V_ull,'U',tk(i).U_N2/tk(i).m_N2);
        tk(i).T_ull = N2('T','D',tk(i).m_N2/tk(i).V_ull,'U',tk(i).U_N2/tk(i).m_N2);
        tk(i).Pv = interpSat(sat{i}, 2, tk(i).T_liq);
        tk(i).P = tk(i).P_N2 + tk(i).Pv; % saturated ullage: N2 + vapor at Pv
        tk(i).rho = rhoFun{i}(tk(i).P, tk(i).T_liq);
    end

    % New COPV state from (density, internal energy)
    P_copv = N2('P','D',m_copv/V_copv,'U',U_copv/m_copv);
    T_copv = N2('T','D',m_copv/V_copv,'U',U_copv/m_copv);
end

idx = 1:k;
t_log = t_log(idx);
P_copv_log = P_copv_log(idx);
T_copv_log = T_copv_log(idx);
m_copv_log = m_copv_log(idx);
Pc_log = Pc_log(idx);
F_log = F_log(idx);
OF_log = OF_log(idx);
OFimp_log = OFimp_log(idx);
TMR_log = TMR_log(idx);

mdot_N2_log = mdot_N2_log(idx,:);
r_log = r_log(idx,:);
P_vout_log = P_vout_log(idx,:);
theta_log = theta_log(idx,:);
thetaCmd_log = thetaCmd_log(idx,:);
thetaFF_log = thetaFF_log(idx,:);
e_log = e_log(idx,:);
Kp_log = Kp_log(idx,:);
Cv_log = Cv_log(idx,:);

P_tank_log = P_tank_log(idx,:);
Pv_log = Pv_log(idx,:);
T_ull_log = T_ull_log(idx,:);
T_liq_log = T_liq_log(idx,:);
m_N2_log = m_N2_log(idx,:);
m_vap_log = m_vap_log(idx,:);
m_liq_log = m_liq_log(idx,:);
V_ull_log = V_ull_log(idx,:);

mdot_liq_log = mdot_liq_log(idx,:);
P_injIn_log = P_injIn_log(idx,:);
dP_inj_log = dP_inj_log(idx,:);

mdot_annulus_log = mdot_liq_log(:,1)*mdot_annulus/mdot_eth_des; % annulus share of the total ethanol flow
stiff_log = dP_inj_log/Pc_des*100; % injector stiffness dP_inj/Pc [%]
superInj_log = P_injIn_log(:,2) - Pv_log(:,2); % nitrous supercharge at the injector inlet [Pa]
Isp_log = F_log./(mdot_liq_log(:,1) + mdot_liq_log(:,2))/g0; % [s]

% Events
inBurn = ~isnan(t_burn) & t_log >= t_burn;
settled = inBurn & t_log >= t_burn + 0.5;
if ~any(settled), settled = inBurn; end % burn shorter than 0.5 s

% Minimum COPV pressure that keeps each valve choked at its burn-start tank pressure
P_copv_chokeMin = (P_set + dP_line)/r_crit; % r_crit at the final COPV state










%% Results
fprintf('=== Setup ===\n');
if eregOn, eregStr = 'ON (regulating through the burn)'; else, eregStr = 'OFF (shut at burn trigger, blowdown burn)'; end
fprintf('E-REG %s: %g Hz, %g ms delay, %g deg/s servo\n', eregStr, f_ctrl, t_delay*1e3, servo_rate);
for i = 1:2
    fprintf('  %-7s Kp = %g + %g/P | theta_FF = %.2f + %.0f/P (%g deg while pressurizing) | limits %g-%g deg | N2 line loss %g psi\n', legName{i}, ...
        kp_0(i), a_kp(i), theta_0(i), b_ff(i), theta_crack, theta_min, theta_max(i), Pa2psi(dP_line(i)));
end
fprintf('  Feedforward fit (isentropic from %g psig, %g F): %s\n', P_copv_0_psig, T_copv_0_F, ...
    strjoin(arrayfun(@(p,T) sprintf('%.0f psig/%.1f K', p, T), P_fit_psig, T_fit_K, 'UniformOutput', false), ', '));
fprintf('Design: Pc = %g psig | ethanol inj inlet %g psig (dP %.0f psi) + feed %g psi -> tank %.1f psig\n', ...
    Pc_psig, P_injIn_eth_psig, Pa2psi(dP_inj_eth_des), dP_ethFeed_psi, psig(P_set_eth));
fprintf('                       | nitrous inj inlet %g psig (dP %.0f psi) + feed %g psi -> tank %.1f psig\n', ...
    P_injIn_ox_psig, Pa2psi(dP_inj_ox_des), dP_oxFeed_psi, psig(P_set_ox));
P_crit_ox = py.CoolProp.CoolProp.PropsSI('Pcrit','NitrousOxide'); % N2O critical pressure [Pa]
if g2a(P_ox_0_psig) < P_crit_ox
    dT_sub_ox = convtemp(NOX('T','P',g2a(P_ox_0_psig),'Q',0),'K','F') - T_ox_0_F; % subcooling [F]
else
    dT_sub_ox = NaN; % above the critical pressure: no saturation temperature, compressed liquid
end
fprintf('Nitrous: %.1f F fill, Pv %.1f psig, start %.1f psig (%.1f psi N2 pad, %.1f F subcooled) -> %.1f psi supercharge at injector inlet\n', ...
    T_ox_0_F, psig(Pv_ox_0), P_ox_0_psig, P_ox_0_psig - psig(Pv_ox_0), dT_sub_ox, P_injIn_ox_psig - psig(Pv_ox_0));
fprintf('         Dyer flow through the pintle, CdA %.2f mm^2 (HEM choke at %.2f Pv)\n', CdA_ox*1e6, r_hem);
fprintf('Pintle injector: nitrous through the pintle (%.4f kg/s, %.1f m/s), ethanol through the annulus (%.4f kg/s, %.1f m/s), film %.4f kg/s\n', ...
    mdot_ox_des, v_pin_des, mdot_annulus, v_ann_des, mdot_film);
fprintf('                 impingement O/F %.2f (overall %.2f incl. film), TMR %.2f\n', OF_imp_des, OF_des, TMR_des);
fprintf('TCA: At %.1f mm^2, c* %.0f m/s (from design point), L* fill lag %.2f ms, Cf %.3f -> design thrust %.0f N (%.0f lbf), Isp %.1f s\n', ...
    At*1e6, cstar_des, tau_c*1e3, F_des/(Pc_des*At), F_des, F_des/4.44822, F_des/mdot_tot_des/g0);
fprintf('Burn trigger: ethanol %.1f psig (setpoint %+g psi), nitrous %.1f psig (setpoint %+g psi)\n', ...
    psig(P_trig(1)), ethTriggerOffset, psig(P_trig(2)), nitTriggerOffset);
for i = 1:2
    fprintf('%s tank: %.2f L total, %.2f L ullage | loaded %.3f kg (%.3f kg residual)\n', legName{i}, ...
        tk(i).V*1e3, tk(i).V_ull_0*1e3, m_0(i), m_resid(i));
end
fprintf('COPV start: %.3f kg N2 | min COPV pressure to stay choked at design tank pressure: eth %.0f psig, ox %.0f psig\n', ...
    m_copv_0, psig(P_copv_chokeMin(1)), psig(P_copv_chokeMin(2)));

fprintf('\n=== Pressurization ===\n');
for i = 1:2
    if isnan(tk(i).t_reach)
        fprintf('%s tank never reached its burn trigger %.1f psig.\n', legName{i}, psig(P_trig(i)));
    else
        fprintf('%s tank reached its burn trigger %.1f psig at t = %.3f s\n', legName{i}, psig(P_trig(i)), tk(i).t_reach);
    end
end
if ~isnan(t_burn), fprintf('Burn started at t = %.3f s\n', t_burn); end
if ~isnan(t_ereg), fprintf('E-REG burn mode at t = %.3f s (%+.3f s from burn start, feedforward ramps over %g s)\n', ...
        t_ereg, t_ereg - (max([tk.t_reach]) + t_hold), t_mav); end

fprintf('\n=== Valve choking ===\n');
for i = 1:2
    iOpen = Cv_log(:,i) > 0; % steps with the valve open
    if any(iOpen)
        fprintf('%s: N2 stayed CHOKED the whole run (max P_valve_out/P_copv %.3f with valve open, r_crit %.3f at end, COPV %.1f psig at end)\n', ...
            legName{i}, max(r_log(iOpen,i)), r_crit, psig(P_copv_log(end)));
    else
        fprintf('%s: valve never opened.\n', legName{i});
    end
end

fprintf('\n=== Burn ===\n');
fprintf('Ended (%s) at t = %.3f s\n', stopReason, t_log(end));
if any(inBurn)
    tb = t_log(inBurn);
    m_used = trapz(tb, mdot_liq_log(inBurn,1) + mdot_liq_log(inBurn,2)); % propellant burned [kg]
    I_tot = trapz(tb, F_log(inBurn)); % total impulse [N-s]
    fprintf('Burn duration: %.3f s (target %g s)\n', t_log(end) - t_burn, burnTime);
    fprintf('Ethanol flow: %.4f - %.4f kg/s (design %.4f), annulus %.4f - %.4f kg/s\n', ...
        min(mdot_liq_log(settled,1)), max(mdot_liq_log(settled,1)), mdot_eth_des, min(mdot_annulus_log(settled)), max(mdot_annulus_log(settled)));
    fprintf('Nitrous flow (pintle): %.4f - %.4f kg/s (design %.4f) | O/F %.2f - %.2f (design %.2f)  [after first 0.5 s]\n', ...
        min(mdot_liq_log(settled,2)), max(mdot_liq_log(settled,2)), mdot_ox_des, min(OF_log(settled)), max(OF_log(settled)), OF_des);
    fprintf('Pintle impingement: O/F %.2f - %.2f (design %.2f) | TMR %.2f - %.2f (design %.2f)  [after first 0.5 s]\n', ...
        min(OFimp_log(settled)), max(OFimp_log(settled)), OF_imp_des, min(TMR_log(settled)), max(TMR_log(settled)), TMR_des);
    fprintf('Chamber: %.1f - %.1f psig | thrust %.0f - %.0f N (%.0f - %.0f lbf) | Isp %.1f - %.1f s  [after first 0.5 s]\n', ...
        psig(min(Pc_log(settled))), psig(max(Pc_log(settled))), min(F_log(settled)), max(F_log(settled)), ...
        min(F_log(settled))/4.44822, max(F_log(settled))/4.44822, min(Isp_log(settled)), max(Isp_log(settled)));
    fprintf('Total impulse %.0f N-s (%.0f lbf-s), average thrust %.0f N, average Isp %.1f s, propellant burned %.3f kg\n', ...
        I_tot, I_tot/4.44822, I_tot/(t_log(end) - t_burn), I_tot/m_used/g0, m_used);
    for i = 1:2
        fprintf('%s tank during burn: %.1f - %.1f psig | injector inlet %.1f - %.1f psig | injector dP %.1f - %.1f psi | stiffness %.0f - %.0f %% of Pc\n', ...
            legName{i}, psig(min(P_tank_log(inBurn,i))), psig(max(P_tank_log(inBurn,i))), ...
            psig(min(P_injIn_log(settled,i))), psig(max(P_injIn_log(settled,i))), ...
            Pa2psi(min(dP_inj_log(settled,i))), Pa2psi(max(dP_inj_log(settled,i))), ...
            min(stiff_log(settled,i)), max(stiff_log(settled,i)));
    end
    fprintf('Nitrous: liquid %.1f -> %.1f K, Pv %.1f -> %.1f psig, min supercharge at injector inlet %.1f psi, %.3f kg evaporated into ullage\n', ...
        T_liq_log(find(inBurn,1),2), T_liq_log(end,2), psig(Pv_log(find(inBurn,1),2)), psig(Pv_log(end,2)), ...
        Pa2psi(min(superInj_log(settled))), m_vap_log(end,2) - m_vap_log(find(inBurn,1),2));
    iB = find(inBurn,1);
    fprintf('N2 used during burn: %.3f kg total (ethanol %.3f kg, nitrous %.3f kg)\n', m_copv_log(iB) - m_copv_log(end), ...
        m_N2_log(end,1) - m_N2_log(iB,1), m_N2_log(end,2) - m_N2_log(iB,2));
    fprintf('\n=== E-REG ===\n');
    for i = 1:2
        [P_peak, iPeak] = max(P_tank_log(:,i));
        fprintf('%s: peak tank pressure %.1f psig at t = %.3f s (setpoint %.1f psig, overshoot %.1f psi)\n', ...
            legName{i}, psig(P_peak), t_log(iPeak), psig(P_set(i)), psig(P_peak) - psig(P_set(i)));
        fprintf('  First 0.5 s of burn: tank %.1f - %.1f psig\n', psig(min(P_tank_log(inBurn & ~settled,i))), psig(max(P_tank_log(inBurn & ~settled,i))));
        if any(settled)
            fprintf('  After 0.5 s: tank %.1f - %.1f psig (mean %.1f), valve %.1f - %.1f deg\n', psig(min(P_tank_log(settled,i))), ...
                psig(max(P_tank_log(settled,i))), psig(mean(P_tank_log(settled,i))), min(theta_log(settled,i)), max(theta_log(settled,i)));
        end
    end
end
fprintf('\n=== Final state ===\n');
fprintf('  COPV:    %.1f psig, %.1f K, %.3f kg N2 (%.0f%% left)\n', psig(P_copv_log(end)), T_copv_log(end), m_copv_log(end), 100*m_copv_log(end)/m_copv_0);
for i = 1:2
    fprintf('  %-7s  ullage %.1f psig, %.1f K, %.3f kg N2, %.2f L | liquid %.2f K, %.3f kg left\n', legName{i}, ...
        psig(P_tank_log(end,i)), T_ull_log(end,i), m_N2_log(end,i), V_ull_log(end,i)*1e3, T_liq_log(end,i), m_liq_log(end,i));
end
fprintf('  Peak tank pressure: ethanol %.1f psig, nitrous %.1f psig\n', ...
    psig(max(P_tank_log(:,1))), psig(max(P_tank_log(:,2))));


%% Plots 
c1 = [0 0.447 0.741];
c2 = [0.85 0.325 0.098];
c3 = [0.466 0.674 0.188];
c4 = [0.494 0.184 0.556];
c5 = [0.301 0.745 0.933];
cLeg = {c2, c5}; % ethanol, nitrous

% --- Figure 1: COPV & N2 supply ---
figure('Name','COPV & N2 supply','Position',[50 60 900 900]);
tl = tiledlayout(4,1,'TileSpacing','compact'); title(tl, 'COPV & N2 Supply');
ax1(1) = nexttile; hold on;
plot(t_log, psig(P_copv_log), 'Color', c1, 'LineWidth', 1.5);
yline(psig(max(P_copv_chokeMin)), ':', 'Min for choked flow', 'Color', [0.4 0.4 0.4]);
ylabel('COPV pressure [psig]');
ax1(2) = nexttile;
plot(t_log, T_copv_log, 'Color', c1, 'LineWidth', 1.5); ylabel('COPV gas temp [K]');
ax1(3) = nexttile;
plot(t_log, m_copv_log, 'Color', c1, 'LineWidth', 1.5); ylabel('COPV N2 mass [kg]');
ax1(4) = nexttile; hold on;
plot(t_log, mdot_N2_log(:,1), 'Color', c2, 'LineWidth', 1.5);
plot(t_log, mdot_N2_log(:,2), 'Color', c5, 'LineWidth', 1.5); ylabel('N2 mdot [kg/s]');
yyaxis right; hold on;
plot(t_log, r_log(:,1), '--', 'Color', c2); plot(t_log, r_log(:,2), '--', 'Color', c5);
ylabel('P_{valve,out}/P_{COPV} (absolute)'); ylim([0 1]); ax1(4).YAxis(2).Color = [0.5 0.5 0.5];
yyaxis left; 
xlabel('Time [s]'); 
legend('Ethanol','Nitrous','Location','best');

% --- Figures 2 & 3: Nitrous and Ethanol tanks ---
axT = {gobjects(0), gobjects(0)};
for i = 1:2
    figure('Name',[legName{i} ' tank'],'Position',[50+50*i 60 900 900]);
    tl = tiledlayout(4,1,'TileSpacing','compact'); title(tl, [legName{i} ' Tank']);
    axT{i}(1) = nexttile; hold on;
    plot(t_log, psig(P_tank_log(:,i)), 'Color', cLeg{i}, 'LineWidth', 1.5);
    plot(t_log, psig(Pv_log(:,i)), ':', 'Color', c4, 'LineWidth', 1.5);
    yline(psig(P_set(i)), ':', 'Target', 'Color', c3);
    ylabel('Pressure [psig]'); legend('Tank','Vapor pressure','Location','best');
    axT{i}(2) = nexttile; hold on;
    plot(t_log, T_ull_log(:,i), 'Color', cLeg{i}, 'LineWidth', 1.5);
    plot(t_log, T_liq_log(:,i), 'Color', c3, 'LineWidth', 1.5);
    ylabel('Temperature [K]'); legend('Ullage N2','Liquid','Location','best');
    axT{i}(3) = nexttile; hold on;
    plot(t_log, m_liq_log(:,i), 'Color', c3, 'LineWidth', 1.5);
    yline(m_resid(i), ':', 'Residual', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
    ylabel('Liquid mass [kg]');
    yyaxis right; hold on;
    plot(t_log, m_N2_log(:,i), '-', 'Color', cLeg{i}, 'LineWidth', 1.5);
    plot(t_log, m_vap_log(:,i), '--', 'Color', c4, 'LineWidth', 1.5);
    ylabel('Ullage mass [kg]'); axT{i}(3).YAxis(2).Color = cLeg{i};
    yyaxis left; legend('Liquid','Ullage N2','Ullage vapor','Location','best');
    axT{i}(4) = nexttile;
    plot(t_log, V_ull_log(:,i)*1e3, 'Color', cLeg{i}, 'LineWidth', 1.5);
    ylabel('Ullage volume [L]'); xlabel('Time [s]');
end

% --- Figure 4: Injector & chamber ---
figure('Name','Injector & chamber','Position',[200 60 900 900]);
tl = tiledlayout(4,1,'TileSpacing','compact'); title(tl, 'Injectors & Chamber');
ax4(1) = nexttile; hold on;
plot(t_log, mdot_liq_log(:,1), 'Color', c2, 'LineWidth', 1.5);
plot(t_log, mdot_annulus_log, '--', 'Color', c2, 'LineWidth', 1.5);
plot(t_log, mdot_liq_log(:,2), 'Color', c5, 'LineWidth', 1.5);
yline([mdot_eth_des mdot_annulus mdot_ox_des], ':', {'Ethanol','Annulus','Nitrous (pintle)'}, 'Color', [0.4 0.4 0.4]);
ylabel('mdot [kg/s]'); legend('Ethanol (tank out)','Ethanol annulus','Nitrous pintle','Location','best');
ax4(2) = nexttile; hold on;
plot(t_log, OF_log, 'Color', c4, 'LineWidth', 1.5);
plot(t_log, OFimp_log, '--', 'Color', c4, 'LineWidth', 1.5);
yline([OF_des OF_imp_des], ':', {'Design','Design (impingement)'}, 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
ylabel('O/F [-]'); ylim([0 2*OF_imp_des]);
yyaxis right; hold on;
plot(t_log, TMR_log, 'Color', c3, 'LineWidth', 1.5);
yline(TMR_des, ':', 'Design TMR', 'Color', c3, 'LabelHorizontalAlignment', 'left', 'HandleVisibility', 'off');
ylabel('TMR [-]'); ax4(2).YAxis(2).Color = c3;
yyaxis left; legend('Overall','Impingement (annulus only)','TMR','Location','best');
ax4(3) = nexttile; hold on;
plot(t_log, psig(P_tank_log(:,1)), 'Color', c2, 'LineWidth', 1);
plot(t_log, psig(P_injIn_log(:,1)), '--', 'Color', c2, 'LineWidth', 1.5);
plot(t_log, psig(P_tank_log(:,2)), 'Color', c5, 'LineWidth', 1);
plot(t_log, psig(P_injIn_log(:,2)), '--', 'Color', c5, 'LineWidth', 1.5);
plot(t_log, psig(Pc_log), 'w', 'LineWidth', 1.5);
yline(Pc_psig, ':', 'Design Pc', 'Color', [0.4 0.4 0.4]);
yline(P_injIn_eth_psig, ':', 'Ethanol inj target', 'Color', c2, 'HandleVisibility', 'off');
yline(P_injIn_ox_psig, ':', 'Nitrous inj target', 'Color', c5, 'HandleVisibility', 'off');
ylabel('Pressure [psig]'); legend('Ethanol tank','Ethanol annulus inlet','Nitrous tank','Nitrous pintle inlet','Chamber','Location','best');
ax4(4) = nexttile; hold on;
plot(t_log, Pa2psi(superInj_log), 'Color', c5, 'LineWidth', 1.5);
yline(0, '--', 'Flashing in injector inlet feed', 'Color', 'r');
ylabel('N2O supercharge at inj [psi]'); xlabel('Time [s]');

% --- Figure 5: Performance ---
figure('Name','Performance','Position',[250 60 900 700]);
tl = tiledlayout(3,1,'TileSpacing','compact'); title(tl, 'Performance');
ax5(1) = nexttile; hold on;
plot(t_log, F_log, 'Color', c1, 'LineWidth', 1.5);
yline(F_des, ':', 'Design', 'Color', [0.4 0.4 0.4]);
ylabel('Thrust [N]');
yyaxis right; ylabel('Thrust [lbf]'); ax5(1).YAxis(2).Color = c1;
ax5(1).YAxis(2).Limits = ax5(1).YAxis(1).Limits/4.44822;
yyaxis left;
ax5(2) = nexttile;
plot(t_log, Isp_log, 'Color', c3, 'LineWidth', 1.5); ylabel('Isp [s]');
ax5(3) = nexttile;
plot(t_log, psig(Pc_log), 'w', 'LineWidth', 1.5); ylabel('Chamber pressure [psig]'); xlabel('Time [s]');

% --- Figure 6: E-REG controllers ---
figure('Name','E-REG controllers','Position',[300 60 1400 900]);
tl = tiledlayout(4,2,'TileSpacing','compact','TileIndexing','columnmajor'); title(tl, 'E-REG Controllers');
ax6 = gobjects(0);
for i = 1:2
    ax6(end+1) = nexttile; hold on;
    plot(t_log, psig(P_tank_log(:,i)), 'Color', cLeg{i}, 'LineWidth', 1.5);
    yline(psig(P_set(i)), '--', 'Setpoint', 'Color', 'w');
    if dP_trig(i) ~= 0, yline(psig(P_trig(i)), ':', 'Burn trigger', 'Color', c3); end
    ylabel([legName{i} ' PT [psig]']); title(legName{i});
    ax6(end+1) = nexttile;
    plot(t_log, e_log(:,i), 'Color', c4, 'LineWidth', 1.5); ylabel('Error e [psi]');
    ax6(end+1) = nexttile; hold on;
    plot(t_log, thetaCmd_log(:,i), ':', 'Color', 'w', 'LineWidth', 1);
    plot(t_log, theta_log(:,i), 'Color', c1, 'LineWidth', 1.5);
    plot(t_log, thetaFF_log(:,i), '--', 'Color', c3, 'LineWidth', 1);
    yline(theta_max(i), ':', 'Limit', 'Color', [0.4 0.4 0.4]);
    ylabel('Valve angle [deg]'); ylim([0 90]); legend('Commanded','Actual','Feedforward','Location','best');
    ax6(end+1) = nexttile; hold on;
    plot(t_log, Cv_log(:,i), 'Color', c1, 'LineWidth', 1.5); ylabel('Valve Cv');
    yyaxis right; plot(t_log, Kp_log(:,i), 'Color', c4); ylabel('K_p(P) [deg/psi]'); ax6(end).YAxis(2).Color = c4;
    yyaxis left; xlabel('Time [s]');
end

% Mark burn start on every plot
markEvents([ax1 axT{1} axT{2} ax4 ax5 ax6], t_log(end), t_burn, t_ereg);

%% Local functions---------------------------------------------------------
function markEvents(axs, t_end, t_burn, t_ereg)
    for a = axs
        axes(a); grid on; xlim([0 t_end]);
        if ~isnan(t_burn)
            xline(t_burn, 'g-', 'Burn start', 'LineWidth', 1.5, 'HandleVisibility', 'off', 'LabelVerticalAlignment', 'bottom');
        end
        if ~isnan(t_ereg) && abs(t_ereg - t_burn) > 1e-6 % only when offset from burn start
            xline(t_ereg, 'm--', 'E-REG burn mode', 'LineWidth', 1.2, 'HandleVisibility', 'off', 'LabelVerticalAlignment', 'bottom');
        end
    end
    linkaxes(axs, 'x');
end

% Saturation table for the propellant vapor in the ullage
% Columns: T [K], Pv [Pa], rho_vap [kg/m^3], h_fg [J/kg], cp_liq [J/kg-K]
function tab = satTable(fluid, T)
    tab = zeros(numel(T), 5);
    for j = 1:numel(T)
        tab(j,:) = [T(j), fluid('P','T',T(j),'Q',0), fluid('D','T',T(j),'Q',1), ...
            fluid('H','T',T(j),'Q',1) - fluid('H','T',T(j),'Q',0), fluid('Cpmass','T',T(j),'Q',0)];
    end
end

% Interpolate column col of a saturation table at T (clamped to the table)
function y = interpSat(tab, col, T)
    y = interp1(tab(:,1), tab(:,col), min(max(T, tab(1,1)), tab(end,1)));
end

% Saturation table for the ethanol/water fuel vapor in the ullage
% Columns: T [K], Pv [Pa], rho_vap [kg/m^3], h_fg [J/kg], cp_liq [J/kg-K]
function tab = satTableFuel(T, x_E, w_E, M_E, M_W, EtOH, H2O)
Ru = 8.314462618; % [J/mol-K]
tab = zeros(numel(T), 5);
    for j = 1:numel(T)
        p_E = x_E*EtOH('P','T',T(j),'Q',0); % ethanol partial pressure [Pa]
        p_W = (1 - x_E)*H2O('P','T',T(j),'Q',0); % water partial pressure [Pa]
        Pv = p_E + p_W;
        y_E = p_E/Pv; % ethanol mole fraction in the vapor
        M_v = y_E*M_E + (1 - y_E)*M_W; % vapor molar mass [kg/mol]
        w_vE = y_E*M_E/M_v; % ethanol mass fraction in the vapor
        h_fg = w_vE*(EtOH('H','T',T(j),'Q',1) - EtOH('H','T',T(j),'Q',0)) + (1 - w_vE)*(H2O('H','T',T(j),'Q',1) - H2O('H','T',T(j),'Q',0));
        cp_liq = w_E*EtOH('Cpmass','T',T(j),'Q',0) + (1 - w_E)*H2O('Cpmass','T',T(j),'Q',0); % ideal-mixing liquid cp
        tab(j,:) = [T(j), Pv, Pv*M_v/(Ru*T(j)), h_fg, cp_liq];
    end
end

% Nitrous mass flux through the pintle per unit CdA [kg/s-m^2]
% Dyer: (k*SPI + HEM)/(1 + k), k = sqrt((P1 - P2)/(Pv - P2))
% SPI: single-phase incompressible (no flashing)
% HEM: isentropic equilibrium expansion to P2, choked below r_hem*Pv
function G = oxFlux(P1, P2, T1, rho1, Pv, r_hem, NOX)
% used: https://ntrs.nasa.gov/api/citations/20190001326/downloads/20190001326.pdf
    dP = max(P1 - P2, 0);
    G_spi = sqrt(2*rho1*dP);
    if P2 >= Pv || dP == 0
        G = G_spi; % liquid can't flash, Dyer reduces to SPI
        return
    end
    P2h = max(P2, r_hem*Pv); % HEM back pressure, limited to below its critical point
    s1 = NOX('S','P',P1,'T',T1);
    h1 = NOX('H','P',P1,'T',T1);
    h2 = NOX('H','P',P2h,'S',s1);
    rho2 = NOX('D','P',P2h,'S',s1);
    G_hem = rho2*sqrt(2*max(h1 - h2, 0));
    k = sqrt(dP/(Pv - P2));
    G = (k*G_spi + G_hem)/(1 + k);
end