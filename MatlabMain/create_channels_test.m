% This script creates a text file where each line contains
% two comma-separated values (lower_limit, upper_limit)

cd channel_select\

% 1. Define the range and number of channels
min_energy = 1;     % Starting value (lower limit of the first bin)
max_energy = 100;   % Ending value (upper limit of the last bin)
num_channels = 50;  % Total number of channels/bins

bin_edges = logspace(log10(min_energy), log10(max_energy), num_channels + 1);

% The lower limits are all edges except the last one.
lower_limits = bin_edges(1:end-1); 
% The upper limits are all edges except the first one.
upper_limits = bin_edges(2:end); 

filename = 'proton_channels_test.txt';

fileID = fopen(filename, 'w');

if fileID == -1
    error('Could not open file for writing: %s', filename);
end

for i = 1 : num_channels
    lower_lim = lower_limits(i);
    upper_lim = upper_limits(i);
    
    % Write the formatted string to the file (using 4 decimal places for precision)
    fprintf(fileID, '%.4f, %.4f\n', lower_lim, upper_lim);
end

fclose(fileID);
cd ../
