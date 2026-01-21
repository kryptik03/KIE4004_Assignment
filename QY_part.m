%% ------------------------------------------------------------------------
%  TASK 3: RE Integration Analysis (Voltage, Losses, Stability)
%  ------------------------------------------------------------------------
%  Objective: Analyze impact on Voltage, Losses, and Stability (FVSI).

fprintf('\n=== Running Task 3: RE Integration Analysis ===\n');

% --- Configuration ---
sys_id_re = 118;          % Test System (e.g., 33, 69, 118)
re_config.enabled = 1;
re_config.bus_idx = 77;   % Location
re_config.p_mw    = 2.5;  % Capacity (MW)
re_config.q_mvar  = 0;    % Unity Power Factor
re_config.type    = 'PQ';

% 1. Build & Run Base Case (No RE)
mpc_base = build_dist_model(sys_id_re);
results_base = runpf(mpc_base, mpoption('verbose', 0, 'out.all', 0));

% 2. Build & Run RE Case (With RE)
mpc_re = build_dist_model(sys_id_re, re_config);
results_re = runpf(mpc_re, mpoption('verbose', 0, 'out.all', 0));

% --- ANALYSIS: Calculate the 3 Key Metrics ---

% 1. Voltage Profile (Minimum Voltage)
v_min_base = min(results_base.bus(:,8));
v_min_re   = min(results_re.bus(:,8));

% 2. System Losses (Total Active Power Loss)
loss_base = sum(real(get_losses(results_base)));
loss_re   = sum(real(get_losses(results_re)));

% 3. Stability (Fast Voltage Stability Index - FVSI)
% We need the helper function 'calc_max_fvsi' at the bottom of the script.
fvsi_base = calc_max_fvsi(results_base);
fvsi_re   = calc_max_fvsi(results_re);

% --- PRINT RESULTS TABLE ---
fprintf('\n--------------------------------------------------\n');
fprintf(' RESULTS SUMMARY (IEEE %d-Bus)\n', sys_id_re);
fprintf(' Integration: %0.2f MW @ Bus %d\n', re_config.p_mw, re_config.bus_idx);
fprintf('--------------------------------------------------\n');
fprintf('Metric           | Base Case  | With RE    | Improvement\n');
fprintf('-----------------|------------|------------|------------\n');
fprintf('1. Min Voltage   | %0.4f pu   | %0.4f pu   | %+0.4f pu\n', v_min_base, v_min_re, v_min_re - v_min_base);
fprintf('2. Total Losses  | %0.4f MW   | %0.4f MW   | %0.2f%%\n', loss_base, loss_re, (loss_re - loss_base)/loss_base * 100);
fprintf('3. Stability     | %0.4f      | %0.4f      | %s\n', fvsi_base, fvsi_re, "Lower is Better");
fprintf('   (Max FVSI)    |            |            | \n');
fprintf('--------------------------------------------------\n');

% --- VISUALIZATION (3 Subplots) ---
figure('Name', ['Task 3 Analysis - IEEE ' num2str(sys_id_re)], 'Color', 'w', 'Position', [100, 50, 600, 800]);

% Subplot 1: Voltage Profile
subplot(3,1,1);
hold on; grid on; grid minor;
plot(results_base.bus(:,1), results_base.bus(:,8), '--r', 'LineWidth', 1.5, 'DisplayName', 'Base Case');
plot(results_re.bus(:,1), results_re.bus(:,8), '-b', 'LineWidth', 1.5, 'DisplayName', 'With RE');
xline(re_config.bus_idx, '-k', 'DisplayName', 'Injection Point');
yline(0.95, '--k', 'DisplayName', 'Limit (0.95)');
title('1. Voltage Profile Impact');
ylabel('Voltage (p.u.)');
legend('Location', 'best');

% Subplot 2: System Losses
subplot(3,1,2);
b1 = bar([loss_base, loss_re], 'FaceColor', 'flat');
b1.CData(1,:) = [1 0 0]; b1.CData(2,:) = [0 0 1];
xticklabels({'Base Case', 'With RE'});
ylabel('Losses (MW)');
title('2. Active Power Loss Reduction');
text(1, loss_base, sprintf('%.3f MW', loss_base), 'Vert', 'bottom', 'Horiz', 'center');
text(2, loss_re,   sprintf('%.3f MW', loss_re),   'Vert', 'bottom', 'Horiz', 'center', 'FontWeight', 'bold');
ylim([0, max(loss_base, loss_re)*1.3]);

% Subplot 3: Stability (FVSI)
subplot(3,1,3);
b2 = bar([fvsi_base, fvsi_re], 'FaceColor', 'flat');
b2.CData(1,:) = [0.85 0.33 0.10]; % Red (High Risk)
b2.CData(2,:) = [0.47 0.67 0.19]; % Green (More Stable)
xticklabels({'Base Case', 'With RE'});
ylabel('Max FVSI Index');
title('3. Stability Improvement (Lower = More Stable)');
text(1, fvsi_base, sprintf('%.3f', fvsi_base), 'Vert', 'bottom', 'Horiz', 'center');
text(2, fvsi_re,   sprintf('%.3f', fvsi_re),   'Vert', 'bottom', 'Horiz', 'center', 'FontWeight', 'bold');
yline(1.0, '--r', 'Unstable Limit (1.0)', 'LineWidth', 1.5);
ylim([0, 1.2]);

% -------------------------------------------------------------------------
%  HELPER FUNCTION: Calculate Stability Index (FVSI)
%  (PASTE THIS AT THE VERY END OF YOUR FILE, AFTER build_dist_model)
% -------------------------------------------------------------------------
function max_fvsi = calc_max_fvsi(results)
    % Calculates the Fast Voltage Stability Index (FVSI) for all lines.
    % Formula: FVSI = (4 * Z^2 * Q_recv) / (V_send^2 * X_line)
    % Stability Condition: FVSI must be < 1.00. 
    % A value close to 1.00 indicates the line is near voltage collapse.

    mpc = results;
    nl = size(mpc.branch, 1);
    fvsi_vals = zeros(nl, 1);
    
    for k = 1:nl
        f = mpc.branch(k, 1); % From Bus ID
        t = mpc.branch(k, 2); % To Bus ID
        r = mpc.branch(k, 3); % Resistance
        x = mpc.branch(k, 4); % Reactance
        z2 = r^2 + x^2;       % Impedance Squared
        
        % Get Sending End Voltage (Magnitude)
        v_i = results.bus(results.bus(:,1)==f, 8); 
        
        % Get Reactive Power at Receiving End (Q_to)
        % MATPOWER Branch Col 14 is Q_to in MVAR.
        % We convert to p.u. by dividing by BaseMVA (usually 100).
        q_j = abs(results.branch(k, 14)) / mpc.baseMVA; 
        
        % Avoid division by zero
        if abs(x) > 1e-5 && ~isempty(v_i)
            fvsi_vals(k) = (4 * z2 * q_j) / (v_i^2 * abs(x));
        else
            fvsi_vals(k) = 0;
        end
    end
    
    % Return the Maximum FVSI found in the system (The most critical line)
    max_fvsi = max(fvsi_vals);
end

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