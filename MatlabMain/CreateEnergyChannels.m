clc
clear all
cd 'D:\HERT_Drive\MATLAB Main\Result'
particle_type = 'proton'; % Change to 'proton' to run the proton case

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
    DE_threshold = 0.6;
    linear_step = 1.0;
    log_break = 53; 
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

%% 2. Generate Linear Channels
Linlowtemp = ExLow_Limit;
current_idx = 1;

while true
    % Calculate remaining bins for the low-energy section
    bins_left = diffLow_total - current_idx + 1;
    
    % Check prospective log bin width if we switched to log right now
    temp_edges = logspace(log10(Linlowtemp), log10(log_break), bins_left + 1);
    temp_DE = temp_edges(2) - temp_edges(1);
    
    if temp_DE < DE_threshold
        % Prospective log bin is too small; create a linear bin
        Channels(current_idx, 1) = Linlowtemp;
        Linlowtemp = Linlowtemp + linear_step;
        Channels(current_idx, 2) = Linlowtemp;
        current_idx = current_idx + 1;
    else
        % Prospective log bin is large enough; break to start log bins
        break;
    end
end

%% 3. Record Core Log Channels
% Generate the continuous, perfectly log-spaced array for the remaining bins
bins_left = diffLow_total - current_idx + 1;
LogLowEdges = logspace(log10(Linlowtemp), log10(log_break), bins_left + 1);

for j = 1:bins_left
    Channels(current_idx, 1) = LogLowEdges(j);
    Channels(current_idx, 2) = LogLowEdges(j+1);
    current_idx = current_idx + 1;
end

%% 4. Generate High Energy Channels
LogHighEdges = logspace(log10(log_break), log10(Max_Limit), diffHigh_total + 1);

for j = 1:diffHigh_total
    Channels(current_idx, 1) = LogHighEdges(j);
    Channels(current_idx, 2) = LogHighEdges(j+1);
    current_idx = current_idx + 1;
end

% Append the final integral bin
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