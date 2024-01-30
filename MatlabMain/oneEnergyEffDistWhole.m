function [singleMatrix_whole, energy_beam, beam_number, count_back_whole, hits_detectors_whole, count_reject] = oneEnergyEffDistWhole(file_name, energy_channels, back_limit, inputfolder, detector_threshold)
% Author: Yinbo Chen
% Date: 6/15/2021
% Modified by: Skyler Krantz, Will Teague
% Date: Nov. 13, 2023

numDetect = 9;

% Initialize output variables
singleMatrix_whole = zeros(size(energy_channels,1), 1);
hits_log = 0;
hits_detectors_whole = zeros(numDetect, 1);

% Find the first sequence of digits and convert it to a number
beam_number = str2double(regexp(file_name, '\d+', 'match', 'once'));
run_number = str2double(regexp(file_name, '\d+_Run(\d+)', 'tokens', 'once'));

% Change to the input folder
cd(inputfolder);

% Opens the file to be processed
fide = fopen(file_name, 'r');
% Sorts each line in the .txt file into a cell array
data = textscan(fide, '%f %f %f %f %f %f %f %f %f %f', 'Delimiter','','HeaderLines',1);  % Skip header
NumNoEnergy = textscan(fide, '%f', 'Delimiter','','HeaderLines',1);
NumEnergyDeposit = textscan(fide, '%f', 'Delimiter','','HeaderLines',1);
% Closes the file
fclose(fide);

% Get the number of rows in the data file
n = size(data{1, 1},1);

% Extract relevant information from the data
energy_beam = data{1};
Detector_Energy = cell2mat(data(2:end));
NumEnergyDeposit = NumEnergyDeposit{1, 1};
NumNoEnergy = NumNoEnergy{1, 1};

% Check if the sum of energy deposits and no-energy counts match the beam number
if NumNoEnergy + NumEnergyDeposit ~= beam_number
    error('ERROR IN THE DATA FILE!');
end

% Reset counters
count_back_whole = 0;

percentage = 0;
fprintf('Starting Run %d \nPercent Complete: ', run_number)
h = 1;

% Zero out values below detector threshold for the whole configuration
for i = 1:numDetect
    if Detector_Energy(i) < detector_threshold
        Detector_Energy(i) = 0;
    end
end

% Calculate the sum of energy deposits for the whole configuration
WholeSum = sum(Detector_Energy); 

% Read data for each detector
for j = 1:(numDetect-1)
    % Update hits count for each detector
    hits_detectors_whole(j) = nnz(Detector_Energy(:,j));
end
    
% Update singleMatrix for each energy channel
for k = 1:length(energy_channels(:, 1))
    for i = 1:numDetect
        if WholeSum(i) >= energy_channels(k, 1) && WholeSum(i) < energy_channels(k, 2)
            singleMatrix_whole(k) = singleMatrix_whole(k) + 1;
        end
    end
end

% Display summary information after running through all simulations
fprintf('Whole Configuration: \nNumber of back hits= %i\nTotal number of hits = %i\n', count_back_whole,sum(singleMatrix_whole));

% Change back to the original directory
cd ..

end