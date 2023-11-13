function [singleMatrix_whole, energy_beam, beam_number, hits_log, count_back_whole, detector_energy_whole, hits_detectors_whole] = oneEnergyEffDistWhole(file_name, energy_channels, back_limit, inputfolder, detector_threshold)
% Author: Yinbo Chen
% Date: 6/15/2021
% Modified by: Skyler Krantz, Will Teague
% Date: Nov. 13, 2023

numDetect = 9;

% Initialize output variables
singleMatrix_whole = zeros(length(energy_channels), 1);
hits_log = 0;
hits_detectors_whole = zeros(numDetect, 1);
detector_energy_whole = zeros(numDetect, 1);

% Recognize the particle energy level from the file name
token = str2double(split(regexp(file_name, '\d+\.*\d*\_\d+', 'match', 'once'), "_"));
energy_beam = token(1);
beam_number = token(2);

% Change to the input folder
cd(inputfolder);

% Read data from the file
fide = fopen(file_name, 'r');
data = textscan(fide, '%s', 'delimiter', '\n');
fclose(fide);

% Get the number of rows in the data file
[n, ~] = size(data{1, 1});

% Extract relevant information from the data
NumEnergyDeposit = str2double(data{1, 1}{n, 1});
NumNoEnergy = str2double(data{1, 1}{n - 2, 1});

% Check if the sum of energy deposits and no-energy counts match the beam number
if NumNoEnergy + NumEnergyDeposit ~= beam_number
    error('ERROR IN THE DATA FILE!');
end

% Reset counters
count_back_whole = 0;

percentage = 0;
fprintf('Starting %.2f \nPercent Complete:', energy_beam)
h = 1;

% Loop through all simulations in the file
for i = 1:NumEnergyDeposit
    % Print percentage complete
    if i > (percentage * (10^-2) * NumEnergyDeposit)
        percentage = percentage + 10;
        fprintf('%.0f ', percentage)
    end
    
    % Reset values for each iteration
    reject_count = 0;

    % Preallocate variables
    Detector_Energy = zeros(numDetect, 1);

    % Read data for each detector
    for j = 1:numDetect
        Detector_Energy(j) = str2double(data{1, 1}{j + 1 + (i - 1) * (numDetect + 1), 1});
        
        % Update hits count for each detector
        if Detector_Energy(j) > 0
            hits_detectors_whole(j) = hits_detectors_whole(j) + 1;
        end
    end

    % Update total detector energy
    detector_energy_whole = detector_energy_whole + [Detector_Energy];

    % Check if the back detector has energy greater than the threshold
    if Detector_Energy(numDetect) > back_limit
        count_back_whole = count_back_whole + 1;
    end

    % Zero out values below detector threshold for the whole configuration
    for k = 1:numDetect
        if Detector_Energy(k) < detector_threshold
            Detector_Energy(k) = 0;
        end
    end

    % Calculate the sum of energy deposits for the whole configuration
    WholeSum = sum(Detector_Energy);

    % Update singleMatrix for each energy channel
    if WholeSum > 0 && Detector_Energy(1) > detector_threshold && reject_count == 0
        for k = 1:length(energy_channels(:, 1))
            if WholeSum >= energy_channels(k, 1) && WholeSum < energy_channels(k, 2)
                singleMatrix_whole(k) = singleMatrix_whole(k) + 1;
                hits_log(h) = i;
                h = h + 1;
            end
        end
    end
end

% Display summary information after running through all simulations
fprintf('\nEnergy Level:%.2i\n', energy_beam)
fprintf('Whole Configuration: \nNumber of back hits= %i\nTotal number of hits = %i\n', count_back_whole, sum(singleMatrix_whole));

% Change back to the original directory
cd ..

end