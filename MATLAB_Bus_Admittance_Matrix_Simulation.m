% =========================================================================
% KIE4004 - Task 3: Fault Analysis (Symmetrical Components)
% System: IEEE 69-Bus Distribution System
% Update: Added Sequence Voltage Prints (V0, V1, V2)
% =========================================================================
clear; clc;

% =========================================================================
% 1. USER CONFIGURATION & SOURCE PARAMETERS
% =========================================================================
S_sc_max  = 100;   % Short Circuit Capacity of Grid (MVA)
c_max     = 1.1;   % Voltage Factor (IEC 60909)
V_gen     = 1.0;   % Generator Voltage (p.u.)

% Ratios for Impedance Calculation
ratio_RX_pos = 0.1; % Ratio R/X for Positive Sequence
ratio_X0_X1  = 1.0; % Ratio X0/X1
ratio_R0_X0  = 0.1; % Ratio R0/X0

% Simulation Settings
fault_bus = 10;    % Target Fault Bus
V_pre     = 1.0;   % Pre-fault voltage (p.u.)
Z_f       = 0;     % Fault Impedance (Bolted)

% =========================================================================
% 2. SOURCE IMPEDANCE CALCULATION
% =========================================================================
% Note: Ensure get_case69.m is in your directory
data = get_case69();
baseMVA = data.baseMVA; 
baseKV = data.baseKV;

% A. Positive/Negative Sequence Source Impedance
S_sc_pu = S_sc_max / baseMVA;
Z_mag_source = (c_max * V_gen^2) / S_sc_pu;
X1_source = Z_mag_source / sqrt(ratio_RX_pos^2 + 1);
R1_source = ratio_RX_pos * X1_source;
Z1_source = R1_source + 1j*X1_source;
Z2_source = Z1_source; 

% B. Zero Sequence Source Impedance
X0_source = ratio_X0_X1 * X1_source;
R0_source = ratio_R0_X0 * X0_source;
Z0_source = R0_source + 1j*X0_source;

fprintf('--- Source Impedance Calculated ---\n');
fprintf('Z1_source: %.4f + j%.4f p.u.\n', real(Z1_source), imag(Z1_source));
fprintf('Z0_source: %.4f + j%.4f p.u.\n', real(Z0_source), imag(Z0_source));
fprintf('-----------------------------------\n');

% =========================================================================
% 3. BUILD LINE SEQUENCE NETWORKS (Y-Bus)
% =========================================================================
branch = data.branch; 
nb = 69; 
Y1 = zeros(nb, nb); Y2 = zeros(nb, nb); Y0 = zeros(nb, nb); 

for k = 1:size(branch, 1)
    f = branch(k, 1); t = branch(k, 2);
    r_pu = branch(k, 3); x_pu = branch(k, 4);
    
    z1 = r_pu + 1j*x_pu; y1 = 1 / z1;
    z2 = z1; y2 = 1 / z2;
    z0 = 3 * z1; y0 = 1 / z0;
    
    Y1(f,t)=-y1; Y1(t,f)=-y1; Y1(f,f)=Y1(f,f)+y1; Y1(t,t)=Y1(t,t)+y1;
    Y2(f,t)=-y2; Y2(t,f)=-y2; Y2(f,f)=Y2(f,f)+y2; Y2(t,t)=Y2(t,t)+y2;
    Y0(f,t)=-y0; Y0(t,f)=-y0; Y0(f,f)=Y0(f,f)+y0; Y0(t,t)=Y0(t,t)+y0;
end

% =========================================================================
% 4. BUILD Z-BUS (NETWORK IMPEDANCE)
% =========================================================================
get_Zbus_slack = @(Y) [0, zeros(1, nb-1); zeros(nb-1, 1), inv(Y(2:end, 2:end))];
Zbus1_net = get_Zbus_slack(Y1);
Zbus2_net = get_Zbus_slack(Y2);
Zbus0_net = get_Zbus_slack(Y0);

% =========================================================================
% 5. FAULT CALCULATION LOGIC
% =========================================================================
Z1_th = Zbus1_net(fault_bus, fault_bus) + Z1_source;
Z2_th = Zbus2_net(fault_bus, fault_bus) + Z2_source;
Z0_th = Zbus0_net(fault_bus, fault_bus) + Z0_source;

Ibase = (baseMVA * 1e6) / (sqrt(3) * baseKV * 1e3); 
Vbase = baseKV / sqrt(3); % Phase-to-neutral base voltage in kV
a = exp(1j * 2 * pi / 3); % Operator a = 1 < 120

fprintf('\n======================================================\n');
fprintf('FAULT ANALYSIS REPORT: IEEE 33-BUS SYSTEM\n');
fprintf('======================================================\n');
fprintf('Base MVA: %.0f MVA  |  Base kV: %.2f kV\n', baseMVA, baseKV);
fprintf('Fault Location: Bus %d\n', fault_bus);
fprintf('S_sc_max: %.0f MVA | c_max: %.2f\n', S_sc_max, c_max);
fprintf('------------------------------------------------------\n');

% --- A. Single Line-to-Ground (SLG) Fault (Phase A) ---
% Current Calculation
I_seq_slg = V_pre / (Z1_th + Z2_th + Z0_th + 3*Z_f);
If_slg_pu = c_max * 3 * I_seq_slg;       
If_slg_amp = abs(If_slg_pu) * Ibase;

fprintf('1. Single Line-to-Ground (SLG) Fault\n');
fprintf('   Fault Current (If): %.4f p.u.  (%.2f Amps)\n', abs(If_slg_pu), If_slg_amp);

% --- VOLTAGE CALCULATION (SLG) ---
% 1. Determine Sequence Currents (I0 = I1 = I2 = If/3)
I0 = If_slg_pu / 3; 
I1 = If_slg_pu / 3; 
I2 = If_slg_pu / 3;

Ia = I0 + I1 + I2; Ib = I0 + (a^2)*I1 + a*I2; Ic = I0 + a*I1 + (a^2)*I2;

% 2. Calculate Bus Sequence Voltages
E_source = V_pre; % Using pre-fault voltage as Source E
V0_seq = 0 - (I0 * Z0_th);
V1_seq = E_source - (I1 * Z1_th);
V2_seq = 0 - (I2 * Z2_th);

% 3. Convert to Phase Voltages
Va = V0_seq + V1_seq + V2_seq;
Vb = V0_seq + (a^2 * V1_seq) + (a * V2_seq);
Vc = V0_seq + (a * V1_seq) + (a^2 * V2_seq);

fprintf('   Sequence Voltages (Mag < Angle):\n');
fprintf('   V0 = %.4f < %.2f deg p.u., ', abs(V0_seq), rad2deg(angle(V0_seq)));
fprintf('   V1 = %.4f < %.2f deg p.u., ', abs(V1_seq), rad2deg(angle(V1_seq)));
fprintf('   V2 = %.4f < %.2f deg p.u.\n', abs(V2_seq), rad2deg(angle(V2_seq)));

fprintf('   Phase Voltages (Mag < Angle):\n');
fprintf('   Va = %.4f < %.2f deg p.u., ', abs(Va), rad2deg(angle(Va)));
fprintf('   Vb = %.4f < %.2f deg p.u., ', abs(Vb), rad2deg(angle(Vb)));
fprintf('   Vc = %.4f < %.2f deg p.u.\n', abs(Vc), rad2deg(angle(Vc)));
fprintf('\n');

fprintf('   Sequence Currents (Mag < Angle):\n');
fprintf('   I0 = %.4f < %.2f deg p.u., ', abs(I0), rad2deg(angle(I0)));
fprintf('   I1 = %.4f < %.2f deg p.u., ', abs(I1), rad2deg(angle(I1)));
fprintf('   I2 = %.4f < %.2f deg p.u.\n', abs(I2), rad2deg(angle(I2)));

fprintf('   Phase Currents (Mag < Angle):\n');
fprintf('   Ia = %.4f < %.2f deg p.u., ', abs(Ia), rad2deg(angle(Ia)));
fprintf('   Ib = %.4f < %.2f deg p.u., ', abs(Ib), rad2deg(angle(Ib)));
fprintf('   Ic = %.4f < %.2f deg p.u.\n', abs(Ic), rad2deg(angle(Ic)));
fprintf('\n');

% --- B. Line-to-Line (LL) Fault (Phase B-C) ---
I1_ll = V_pre / (Z1_th + Z2_th + Z_f);
I2_ll = -I1_ll;
I0_ll = 0;

Ib_ll_pu = c_max * (I0_ll + (a^2)*I1_ll + a*I2_ll);
Ib_ll_amp = abs(Ib_ll_pu) * Ibase;

fprintf('2. Line-to-Line (LL) Fault\n');
fprintf('   Fault Current (Ib): %.4f p.u.  (%.2f Amps)\n', abs(Ib_ll_pu), Ib_ll_amp);

% --- VOLTAGE CALCULATION (LL) ---
I1 = c_max * I1_ll; 
I2 = c_max * I2_ll; 
I0 = 0;

Ia = I0 + I1 + I2; Ib = I0 + (a^2)*I1 + a*I2; Ic = I0 + a*I1 + (a^2)*I2;

V0_seq = 0 - (I0 * Z0_th);
V1_seq = E_source - (I1 * Z1_th);
V2_seq = 0 - (I2 * Z2_th);

Va = V0_seq + V1_seq + V2_seq;
Vb = V0_seq + (a^2 * V1_seq) + (a * V2_seq);
Vc = V0_seq + (a * V1_seq) + (a^2 * V2_seq);

fprintf('   Sequence Voltages (Mag < Angle):\n');
fprintf('   V0 = %.4f < %.2f deg p.u., ', abs(V0_seq), rad2deg(angle(V0_seq)));
fprintf('   V1 = %.4f < %.2f deg p.u., ', abs(V1_seq), rad2deg(angle(V1_seq)));
fprintf('   V2 = %.4f < %.2f deg p.u.\n', abs(V2_seq), rad2deg(angle(V2_seq)));

fprintf('   Phase Voltages (Mag < Angle):\n');
fprintf('   Va = %.4f < %.2f deg p.u., ', abs(Va), rad2deg(angle(Va)));
fprintf('   Vb = %.4f < %.2f deg p.u., ', abs(Vb), rad2deg(angle(Vb)));
fprintf('   Vc = %.4f < %.2f deg p.u.\n', abs(Vc), rad2deg(angle(Vc)));
fprintf('\n');

fprintf('   Sequence Currents (Mag < Angle):\n');
fprintf('   I0 = %.4f < %.2f deg p.u., ', abs(I0), rad2deg(angle(I0)));
fprintf('   I1 = %.4f < %.2f deg p.u., ', abs(I1), rad2deg(angle(I1)));
fprintf('   I2 = %.4f < %.2f deg p.u.\n', abs(I2), rad2deg(angle(I2)));

fprintf('   Phase Currents (Mag < Angle):\n');
fprintf('   Ia = %.4f < %.2f deg p.u., ', abs(Ia), rad2deg(angle(Ia)));
fprintf('   Ib = %.4f < %.2f deg p.u., ', abs(Ib), rad2deg(angle(Ib)));
fprintf('   Ic = %.4f < %.2f deg p.u.\n', abs(Ic), rad2deg(angle(Ic)));
fprintf('\n');



% --- C. Double Line-to-Ground (DLG) Fault (Phase B-C-G) ---
Z_eq = (Z2_th * (Z0_th + 3*Z_f)) / (Z2_th + Z0_th + 3*Z_f);
I1_dlg_raw = V_pre / (Z1_th + Z_eq);

% Scaling currents by c_max
I1 = c_max * I1_dlg_raw;
V1_f_temp = E_source - I1 * Z1_th; % Voltage at fault point
I2 = -c_max*V1_f_temp / Z2_th;
I0 = -c_max*V1_f_temp / (Z0_th + 3*Z_f);

Ia = I0 + I1 + I2; Ib = I0 + (a^2)*I1 + a*I2; Ic = I0 + a*I1 + (a^2)*I2;

Ig_dlg_pu = 3 * I0;
Ig_dlg_amp = abs(Ig_dlg_pu) * Ibase;

fprintf('3. Double Line-to-Ground (DLG) Fault\n');
fprintf('   Ground Current (Ig): %.4f p.u.  (%.2f Amps)\n', abs(Ig_dlg_pu), Ig_dlg_amp);

% --- VOLTAGE CALCULATION (DLG) ---
V0_seq = 0 - (I0 * Z0_th);
V1_seq = E_source - (I1 * Z1_th);
V2_seq = 0 - (I2 * Z2_th);

Va = V0_seq + V1_seq + V2_seq;
Vb = V0_seq + (a^2 * V1_seq) + (a * V2_seq);
Vc = V0_seq + (a * V1_seq) + (a^2 * V2_seq);

fprintf('   Sequence Voltages (Mag < Angle):\n');
fprintf('   V0 = %.4f < %.2f deg p.u., ', abs(V0_seq), rad2deg(angle(V0_seq)));
fprintf('   V1 = %.4f < %.2f deg p.u., ', abs(V1_seq), rad2deg(angle(V1_seq)));
fprintf('   V2 = %.4f < %.2f deg p.u.\n', abs(V2_seq), rad2deg(angle(V2_seq)));

fprintf('   Phase Voltages (Mag < Angle):\n');
fprintf('   Va = %.4f < %.2f deg p.u., ', abs(Va), rad2deg(angle(Va)));
fprintf('   Vb = %.4f < %.2f deg p.u., ', abs(Vb), rad2deg(angle(Vb)));
fprintf('   Vc = %.4f < %.2f deg p.u.\n', abs(Vc), rad2deg(angle(Vc)));
fprintf('\n');

fprintf('   Sequence Currents (Mag < Angle):\n');
fprintf('   I0 = %.4f < %.2f deg p.u., ', abs(I0), rad2deg(angle(I0)));
fprintf('   I1 = %.4f < %.2f deg p.u., ', abs(I1), rad2deg(angle(I1)));
fprintf('   I2 = %.4f < %.2f deg p.u.\n', abs(I2), rad2deg(angle(I2)));

fprintf('   Phase Currents (Mag < Angle):\n');
fprintf('   Ia = %.4f < %.2f deg p.u., ', abs(Ia), rad2deg(angle(Ia)));
fprintf('   Ib = %.4f < %.2f deg p.u., ', abs(Ib), rad2deg(angle(Ib)));
fprintf('   Ic = %.4f < %.2f deg p.u.\n', abs(Ic), rad2deg(angle(Ic)));
fprintf('\n');
fprintf('======================================================\n');