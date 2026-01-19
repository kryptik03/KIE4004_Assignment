function case33bw = get_case33bw()
% GET_CASE33BW Returns the IEEE 33-bus (Baran & Wu) system data.
%
% Constants:
%   Base MVA: 10 MVA
%   Base kV:  12.66 kV
%   Z_base:   16.0275 Ohms
%
% Returns:
%   case33bw.branch: [From, To, R(pu), X(pu)]
%   case33bw.baseMVA: 10
%   case33bw.baseKV: 12.66

    case33bw.baseMVA = 10;   
    case33bw.baseKV = 12.66; 

    % =====================================================================
    % 1. RAW BRANCH DATA (Values in OHMS)
    % Cols: [From, To, R(ohm), X(ohm)]
    % Source: Baran and Wu (1989)
    % =====================================================================
    branch_ohms = [
        1   2   0.0922  0.0470;
        2   3   0.4930  0.2511;
        3   4   0.3660  0.1864;
        4   5   0.3811  0.1941;
        5   6   0.8190  0.7070;
        6   7   0.1872  0.6188;
        7   8   0.7114  0.2351;
        8   9   1.0300  0.7400;
        9   10  1.0440  0.7400;
        10  11  0.1966  0.0650;
        11  12  0.3744  0.1238;
        12  13  1.4680  1.1550;
        13  14  0.5416  0.7129;
        14  15  0.5910  0.5260;
        15  16  0.7463  0.5450;
        16  17  1.2890  1.7210;
        17  18  0.7320  0.5740;
        2   19  0.1640  0.1565;
        19  20  1.5042  1.3554;
        20  21  0.4095  0.4784;
        21  22  0.7089  0.9373;
        3   23  0.4512  0.3083;
        23  24  0.8980  0.7091;
        24  25  0.8960  0.7011;
        6   26  0.2030  0.1034;
        26  27  0.1938  0.2259;
        27  28  1.4480  1.5563;
        28  29  1.0590  0.9337;
        29  30  0.5075  0.2585;
        30  31  0.3105  0.3619;
        31  32  0.3410  0.5302;
        32  33  0.6743  0.5403;
    ];

    % =====================================================================
    % 2. CONVERSION: OHMS -> PER UNIT
    % =====================================================================
    
    % Calculate Base Impedance
    % Zbase = (12.66^2) / 10 = 16.02756 Ohm
    Zbase = (case33bw.baseKV^2) / case33bw.baseMVA; 
    
    % Perform conversion
    branch_pu = branch_ohms;
    branch_pu(:, 3) = branch_ohms(:, 3) / Zbase; % R conversion
    branch_pu(:, 4) = branch_ohms(:, 4) / Zbase; % X conversion
    
    case33bw.branch = branch_pu;
end