%% ------------------------------------------------------------------------
%  TASK 1 & 2: Base Case Analysis & Solver Comparison
%  ------------------------------------------------------------------------
%  Objective: Load standard IEEE models and compare NR vs FD methods.
clc; clear; close all;

% Select System ID: 33, 69, or 118
sys_id = 69; 

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

% --- 5. Base Case Results & Power Balance ---
% Extract Real Power Data
P_gen_total  = sum(results_nr.gen(:, 2));      % Total Active Generation (MW)
P_load_total = sum(results_nr.bus(:, 3));      % Total Active Load (MW)
P_loss_total = sum(real(get_losses(results_nr))); % Total Active Losses (MW)

% Extract Voltage Data
[min_volt, min_bus_idx] = min(results_nr.bus(:, 8)); % Min Voltage (p.u.)

% CALCULATE BALANCE (Requirement from Task 1)
% Concept: Generation should equal Load + Losses
balance_mismatch = P_gen_total - P_load_total - P_loss_total;

fprintf('\n=== Base Case Results (IEEE %d-Bus) ===\n', sys_id);
fprintf('Total Generation: %10.4f MW\n', P_gen_total);
fprintf('Total Load:       %10.4f MW\n', P_load_total);
fprintf('Total Losses:     %10.4f MW\n', P_loss_total);
fprintf('---------------------------------------\n');
fprintf('Power Balance:    %10.4f = %10.4f + %10.4f\n', P_gen_total, P_load_total, P_loss_total);
fprintf('Net Mismatch:     %10.6f MW  <-- (Accuracy Check)\n', balance_mismatch);
fprintf('---------------------------------------\n');
fprintf('Min Voltage:      %0.4f p.u. @ Bus %d\n', min_volt, results_nr.bus(min_bus_idx, 1));

% 5. Base Case Results (Voltage & Losses)
total_losses_mw = sum(real(get_losses(results_nr)));
min_volt = min(results_nr.bus(:, 8)); % Column 8 is VM
fprintf('\nBase Case Results:\n');
fprintf('Total System Losses: %0.4f MW\n', total_losses_mw);
fprintf('Minimum Bus Voltage: %0.4f p.u.\n', min_volt);

% --- Visualization: Voltage Profile ---
figure('Name', ['Base Case Voltage - IEEE ' num2str(sys_id)], 'Color', 'w');
plot(results_nr.bus(:,1), results_nr.bus(:,8), '-b', 'LineWidth', 1.5, 'DisplayName', 'Base Case Voltage');
hold on;

% Add Grid Limits (Visual Engineering Insight)
yline(0.95, '--r', 'Min Limit (0.95 p.u.)', 'LineWidth', 1.2, 'DisplayName', 'Grid Code Limit');
yline(1.05, '--r', 'HandleVisibility', 'off'); % Upper limit usually not shown if not violated

% Highlight the Minimum Voltage Bus
plot(results_nr.bus(min_bus_idx, 1), min_volt, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Critical Bus');

title(['Base Case Voltage Profile: IEEE ' num2str(sys_id) '-Bus System'], 'FontSize', 12);
xlabel('Bus Number', 'FontSize', 10);
ylabel('Voltage Magnitude (p.u.)', 'FontSize', 10);
legend('Location', 'best');
grid on;
box on;

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