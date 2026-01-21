%% KIE4004 Power System Assignment - Master Script
% Title: Modeling and Analysis of IEEE Test Systems with RE Integration
% 
% INSTRUCTIONS:
% 1. Ensure MATPOWER is installed and added to your path.
% 2. Run specific tasks by clicking inside a section and pressing "Run Section" (Ctrl+Enter).
% 3. Do not edit the helper function at the bottom unless necessary.

%% ------------------------------------------------------------------------
%  TASK 1 & 2: Base Case Analysis & Solver Comparison
%  ------------------------------------------------------------------------
%  Objective: Load standard IEEE models and compare NR vs FD methods.
clc; clear; close all;

% Select System ID: 33, 69, or 118
sys_id = 33; 

fprintf('=== Running Task 1 & 2: Base Case Analysis (IEEE %d-Bus) ===\n', sys_id);

% 1. Build Base Model (No RE)
mpc_base = build_dist_model(sys_id);

% 2. Run Newton-Raphson (NR)
fprintf('\n--- Newton-Raphson (NR) Method ---\n');
opt_nr = mpoption('pf.alg', 'NR', 'verbose', 1, 'out.all', 0);
t_start = tic;
results_nr = runpf(mpc_base, opt_nr);
t_nr = toc(t_start);

% 3. Run Fast Decoupled (FD)
fprintf('\n--- Fast Decoupled (FDXB) Method ---\n');
opt_fd = mpoption('pf.alg', 'FDXB', 'verbose', 1, 'out.all', 0);
t_start = tic;
results_fd = runpf(mpc_base, opt_fd);
t_fd = toc(t_start);

% 4. Compare Performance
fprintf('\n--- Performance Comparison ---\n');
fprintf('Method | Iterations | Time (sec) | Converged?\n');
fprintf('NR     | %10d | %10.4f | %d\n', results_nr.iterations, t_nr, results_nr.success);
fprintf('FD     | %10d | %10.4f | %d\n', results_fd.iterations, t_fd, results_fd.success);

% 5. Base Case Results (Voltage & Losses)
total_losses_mw = sum(real(get_losses(results_nr)));
min_volt = min(results_nr.bus(:, 8)); % Column 8 is VM
fprintf('\nBase Case Results:\n');
fprintf('Total System Losses: %0.4f MW\n', total_losses_mw);
fprintf('Minimum Bus Voltage: %0.4f p.u.\n', min_volt);


%% ------------------------------------------------------------------------
%  TASK 3: Renewable Energy (RE) Integration
%  ------------------------------------------------------------------------
%  Objective: Inject RE at a specific bus and analyze impact.
%  You can change 're_config' values here to test different scenarios.

fprintf('\n=== Running Task 3: RE Integration Analysis ===\n');

% --- Configuration ---
sys_id_re = 33;          % Test System (e.g., 33)
re_config.enabled = 1;
re_config.bus_idx = 18;  % <--- Change this to test locations
re_config.p_mw    = 2.5; % <--- Change this to test capacity (MW)
re_config.q_mvar  = 0;   % Unity Power Factor
re_config.type    = 'PQ';% Modeling type

% 1. Build Modified Model
mpc_re = build_dist_model(sys_id_re, re_config);

% 2. Run Power Flow
results_re = runpf(mpc_re, mpoption('verbose', 0, 'out.all', 0));

% 3. Compare with Base Case
mpc_base_re = build_dist_model(sys_id_re);
results_base_re = runpf(mpc_base_re, mpoption('verbose', 0, 'out.all', 0));

loss_base = sum(real(get_losses(results_base_re)));
loss_re   = sum(real(get_losses(results_re)));
v_base    = results_base_re.bus(re_config.bus_idx, 8);
v_re      = results_re.bus(re_config.bus_idx, 8);

fprintf('Integration at Bus: %d | Capacity: %0.2f MW\n', re_config.bus_idx, re_config.p_mw);
fprintf('Voltage Impact: %0.4f -> %0.4f p.u.\n', v_base, v_re);
fprintf('Loss Impact:    %0.4f -> %0.4f MW (%0.2f%% Change)\n', ...
    loss_base, loss_re, (loss_re - loss_base)/loss_base * 100);

% Optional: Plot Voltage Profile
figure;
plot(results_base_re.bus(:,8), '-o', 'LineWidth', 1.5, 'DisplayName', 'Base Case');
hold on;
plot(results_re.bus(:,8), '-^', 'LineWidth', 1.5, 'DisplayName', 'With RE');
title(['Voltage Profile (IEEE ' num2str(sys_id_re) '-Bus)']);
xlabel('Bus Index'); ylabel('Voltage (p.u.)');
legend; grid on;


%% ------------------------------------------------------------------------
%  TASK 4 & 5: Fault Analysis Prep (Data Extraction)
%  ------------------------------------------------------------------------
%  Objective: Extract Z matrices and load data for fault calculation.

fprintf('\n=== Running Task 4/5: Data Extraction for Fault Analysis ===\n');

% 1. Load System
sys_id_fault = 69; 
mpc_fault = build_dist_model(sys_id_fault);

% 2. Extract Branch Impedances (Z1 = R + jX)
% Note: Positive sequence is assumed equal to Negative sequence Z1=Z2
% Zero sequence (Z0) usually requires assumption (e.g., Z0 = 3*Z1) for distribution lines.
branches = mpc_fault.branch;
R = branches(:, 3); % Resistance
X = branches(:, 4); % Reactance
Z_lines = R + 1j*X;

% 3. Display Data for Manual Calculation
fprintf('First 5 lines of IEEE %d-Bus System Impedance (Z1):\n', sys_id_fault);
for k = 1:5
    fprintf('Line %d-%d: Z = %0.4f + j%0.4f\n', ...
        branches(k, 1), branches(k, 2), real(Z_lines(k)), imag(Z_lines(k)));
end
fprintf('... (Use "mpc_fault.branch" to see full data)\n');


%% ------------------------------------------------------------------------
%  HELPER FUNCTION: Modular Model Builder
%  ------------------------------------------------------------------------
%  (Do not place any code below this function definition)

function mpc = build_dist_model(sys_id, re_config)
    % BUILD_DIST_MODEL Loads and modifies IEEE distribution systems.
    % Inputs: sys_id (33, 69, 118), re_config (struct)
    
    % 1. Load the Base Case
    switch sys_id
        case 33
            mpc = loadcase('case33bw'); 
        case 69
            mpc = loadcase('case69');   
        case 118
            mpc = loadcase('case118zh'); 
        otherwise
            error('Invalid System ID. Choose 33, 69, or 118.');
    end
    
    % 2. Modular RE Integration
    if nargin > 1 && isfield(re_config, 'enabled') && re_config.enabled
        
        % Constants for MATPOWER indexing
        GEN_BUS = 1; PG = 2; QG = 3; QMAX = 4; QMIN = 5; 
        VG = 6; MBASE = 7; GEN_STATUS = 8; PMAX = 9; PMIN = 10;
        BUS_TYPE = 2;
        
        % --- A. Add New Generator Row ---
        new_gen = zeros(1, 21); 
        new_gen(GEN_BUS)    = re_config.bus_idx;
        new_gen(PG)         = re_config.p_mw;
        new_gen(QG)         = re_config.q_mvar;
        new_gen(GEN_STATUS) = 1;
        new_gen(MBASE)      = mpc.baseMVA;
        
        % Set Limits
        new_gen(PMAX)       = re_config.p_mw;     
        new_gen(PMIN)       = 0;
        new_gen(QMAX)       = re_config.q_mvar + 10;
        new_gen(QMIN)       = re_config.q_mvar - 10;
        new_gen(VG)         = 1.0; % Default voltage setpoint
        
        % Append to gen matrix
        mpc.gen = [mpc.gen; new_gen];
        
        % --- B. Add Corresponding Cost Row (CRITICAL FIX) ---
        % If we don't add this, ext2int crashes because gen and gencost sizes mismatch.
        if isfield(mpc, 'gencost')
            % Create a zero-cost polynomial row: 
            % [Model(2=poly), Startup(0), Shutdown(0), N(3 params), c2(0), c1(0), c0(0)]
            new_gencost = [2, 0, 0, 3, 0, 0, 0];
            mpc.gencost = [mpc.gencost; new_gencost];
        end

        % --- C. Update Bus Type ---
        % If type is 'PV', we must change the bus type to 2 (PV)
        % For 'PQ' (negative load), keeping it as is (usually 1=PQ) is fine, 
        % but standard GEN usually implies PV capability unless constrained.
        if isfield(re_config, 'type') && strcmp(re_config.type, 'PV')
            mpc.bus(re_config.bus_idx, BUS_TYPE) = 2; 
        end
    end
end