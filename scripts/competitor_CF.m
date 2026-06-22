function comp = competitor_CF(atm)
% COMPETITOR_CF  Capacity Factor and AEP for competitor turbines at a given site.
%
%   comp = competitor_CF(atm)
%
%   Input:
%     atm  - struct from wind_resource(), containing at least:
%              .vw_range, .v_mean, .h_ref, .alpha, .k
%
%   Output:
%     comp - struct with one sub-struct per turbine, each containing:
%              .CF   [-]       capacity factor
%              .AEP  [kWh/yr]  annual energy production
%              .P_r  [W]       nameplate power
%              .name           turbine label for display

vw = atm.vw_range;
dv = vw(2) - vw(1);

%% Helper: Weibull PDF at a given hub height
    function gw = weibull_at_h(h_hub)
        v_hub = atm.v_mean * (h_hub / atm.h_ref)^atm.alpha;
        A_hub = v_hub / gamma(1 + 1/atm.k);
        gw    = atm.k/A_hub * (vw/A_hub).^(atm.k-1) .* exp(-(vw/A_hub).^atm.k);
    end

%% Helper: CF from a normalised power curve and hub height
    function CF = compute_CF(P_norm, h_hub)
        gw = weibull_at_h(h_hub);
        CF = sum(P_norm .* gw) * dv;
    end

%% ---- KiteX TWT-11 -------------------------------------------------------
% Source: KiteX official specs + physical estimate
% Rotor Ø11m (A=95m²), cut-in 3.5 m/s, cut-out 18 m/s, P_r=5kW
% Rated speed ~7.0 m/s estimated from Cp≈0.30
h_hub_K  = 15.5;   % [m] official tower height
v_ci_K   = 3.5;    % [m/s] cut-in  (official)
v_r_K    = 7.0;    % [m/s] rated   (estimated)
v_co_K   = 18;     % [m/s] cut-out (official)
P_r_K    = 5e3;    % [W]

P_norm_K = zeros(size(vw));
P_norm_K(vw >= v_ci_K & vw < v_r_K)  = (vw(vw >= v_ci_K & vw < v_r_K).^3 - v_ci_K^3) / (v_r_K^3 - v_ci_K^3);
P_norm_K(vw >= v_r_K  & vw <= v_co_K) = 1;

comp.KiteX.name = 'KiteX TWT-11';
comp.KiteX.P_r  = P_r_K;
comp.KiteX.CF   = compute_CF(P_norm_K, h_hub_K);
comp.KiteX.AEP  = comp.KiteX.CF * P_r_K * 8760 / 1e3;   % [kWh/yr]

%% ---- Bergey Excel 10 ----------------------------------------------------
% Source: SWCC/USDA official test, IEC 61400-12 (Bushland TX, rho=1.225 kg/m³)
% AutoFurl from ~13 m/s; peak 12.6 kW at ~16.5 m/s; no conventional cut-out
% AWEA rated 8.9 kW at 11 m/s; nameplate 10 kW
h_hub_E  = 30;     % [m] tower height (official)
P_r_E    = 10e3;   % [W] nameplate

v_tab_E = [ 0,  1,  2,  3,    4,    5,    6,    7,    8,    9, ...
           10,  11,   12,   13,   14,   15,   16,   17,   18,   19,   20,  21]';
P_tab_E = [ 0,  0,  0,  0,  200,  500, 1600, 2700, 3700, 5700, ...
         7700,10800,12000,12400,12500,12500,12500,12500,12500,12000,11500,11500]';

P_norm_E = max(0, interp1(v_tab_E, P_tab_E/P_r_E, vw, 'pchip', 0));

comp.Excel10.name = 'Bergey Excel 10';
comp.Excel10.P_r  = P_r_E;
comp.Excel10.CF   = compute_CF(P_norm_E, h_hub_E);
comp.Excel10.AEP  = comp.Excel10.CF * P_r_E * 8760 / 1e3;  % [kWh/yr]

%% ---- Bergey Excel 15 ----------------------------------------------------
% Source: official Bergey power curve graph (sea-level, rho=1.225 kg/m³)
% Rotor Ø9.6m; AWEA rated 15.6 kW at 11 m/s; peak ~21 kW at 17 m/s
% AutoFurl above ~15 m/s; no conventional cut-out
h_hub_E15 = 30;     % [m] tower height (standard Bergey tower, same as Excel 10)
P_r_E15   = 15e3;   % [W] nameplate (AWEA rated 15.6 kW)

v_tab_E15 = [ 0,  1,  2,    3,    4,    5,    6,    7,    8,    9, ...
             10,   11,   12,   13,   14,   15,   16,   17,   18,   19]';
P_tab_E15 = [ 0,  0,  0,  200, 1000, 2500, 4500, 7000, 9500,12000, ...
           14500,15600,16500,18000,19000,20000,20500,21000,20000,20000]';

P_norm_E15 = max(0, interp1(v_tab_E15, P_tab_E15/P_r_E15, vw, 'pchip', 0));

comp.Excel15.name = 'Bergey Excel 15';
comp.Excel15.P_r  = P_r_E15;
comp.Excel15.CF   = compute_CF(P_norm_E15, h_hub_E15);
comp.Excel15.AEP  = comp.Excel15.CF * P_r_E15 * 8760 / 1e3;  % [kWh/yr]

end
