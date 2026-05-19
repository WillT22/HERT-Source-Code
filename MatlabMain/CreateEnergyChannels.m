clc
clear all
cd 'D:\HERT_Drive\MATLAB Main\Result'
particle_type = 'proton'; % Change to 'electron' to run the electron case

%% Parameters
if strcmp(particle_type, 'electron')
    ExLow_Limit = 0.1;
    DE_threshold = 0.025;
    linear_step = 0.2;
    log_break = 2.5; 
    Max_Limit = 7.0;
    diffLow_total = 30;
    diffHigh_total = 10;
    
elseif strcmp(particle_type, 'proton')
    ExLow_Limit = 0.1;
    DE_threshold = 0.9;
    linear_step = 0.9;
    % log_break is removed; continuous log interpolation targets Max_Limit
    Max_Limit = 65;
    diffLow_total = 25;
    diffHigh_total = 5;
    
else
    error('Invalid particle type. Choose ''electron'' or ''proton''.');
end

%% 1. Initialize Channels
% Total channels includes low, high, and 1 integral bin at the end
NumChannels = diffLow_total + diffHigh_total + 1;
Channels = zeros(NumChannels, 2);
current_idx = 1;
Linlowtemp = ExLow_Limit;

%% 2-4. Generate Channels (Branched Logic)
if strcmp(particle_type, 'electron')
    % Bound the loop by the total low-energy bins allowed
    while current_idx <= diffLow_total
        bins_left = diffLow_total - current_idx + 1;
        temp_edges = logspace(log10(Linlowtemp), log10(log_break), bins_left + 1);
        temp_DE = temp_edges(2) - temp_edges(1);
        
        if temp_DE < DE_threshold
            Channels(current_idx, 1) = Linlowtemp;
            Linlowtemp = Linlowtemp + linear_step;
            Channels(current_idx, 2) = Linlowtemp;
            current_idx = current_idx + 1;
        else
            break;
        end
    end

    % Record Core Log Channels
    bins_left = diffLow_total - current_idx + 1;
    if bins_left > 0
        LogLowEdges = logspace(log10(Linlowtemp), log10(log_break), bins_left + 1);
        for j = 1:bins_left
            Channels(current_idx, 1) = LogLowEdges(j);
            Channels(current_idx, 2) = LogLowEdges(j+1);
            current_idx = current_idx + 1;
        end
    end

    % Generate High Energy Channels
    LogHighEdges = logspace(log10(log_break), log10(Max_Limit), diffHigh_total + 1);
    for j = 1:diffHigh_total
        Channels(current_idx, 1) = LogHighEdges(j);
        Channels(current_idx, 2) = LogHighEdges(j+1);
        current_idx = current_idx + 1;
    end

elseif strcmp(particle_type, 'proton')
    total_diff_bins = diffLow_total + diffHigh_total;
    
    % Bound the loop by the total differential bins allowed
    while current_idx <= total_diff_bins
        bins_left = total_diff_bins - current_idx + 1;
        
        % Check prospective log bin width targeting Max_Limit directly
        temp_edges = logspace(log10(Linlowtemp), log10(Max_Limit), bins_left + 1);
        temp_DE = temp_edges(2) - temp_edges(1);
        
        if temp_DE < DE_threshold
            Channels(current_idx, 1) = Linlowtemp;
            Linlowtemp = Linlowtemp + linear_step;
            Channels(current_idx, 2) = Linlowtemp;
            current_idx = current_idx + 1;
        else
            break; % Ready to switch to continuous log spacing
        end
    end
    
    % Record Single Continuous Log Array
    bins_left = total_diff_bins - current_idx + 1;
    if bins_left > 0
        LogEdges = logspace(log10(Linlowtemp), log10(Max_Limit), bins_left + 1);
        
        for j = 1:bins_left
            Channels(current_idx, 1) = LogEdges(j);
            Channels(current_idx, 2) = LogEdges(j+1);
            current_idx = current_idx + 1;
        end
    end
end

% Append the final integral bin (Shared Logic)
Channels(current_idx, 1) = Max_Limit;
Channels(current_idx, 2) = 1000;

%% 5. Calculate Energy Resolution
% Calculate resolution for all except the integral bin
Resolution = zeros(NumChannels - 1, 1);
for i = 1:(NumChannels - 1)
   Resolution(i, 1) = 100 * (Channels(i,2) - Channels(i,1)) / ((Channels(i,1) + Channels(i,2)) / 2);
end

disp('Energy Resolutions (%):')
disp(Resolution)

%% 6. Write to Text File
% Dynamically name the output file
filename = sprintf('channel_select\\%s_channels_v6.txt', particle_type);
fileID = fopen(filename, 'w');
if fileID == -1
    error('Could not open file. Ensure the ''channel_select'' directory exists.');
end

for i = 1:size(Channels, 1)
    fprintf(fileID, '%6.3f,%6.3f \n', Channels(i,1), Channels(i,2));
end
fclose(fileID);
fprintf('Data successfully written to %s\n', filename);