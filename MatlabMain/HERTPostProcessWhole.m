%% %% Authored by Skyler Krantz
% Edited by Will Teague 
% Updated: Nov. 13th, 2023
% Post Processing for GEANT4

function [index_arr] = HERTPostProcessWhole(file_name, inputfolder, outputfolder)
% HERTPostProcess Removes simulations with no energy deposited
numDetect = 9;

cd(inputfolder);

% Recognize the particle energy level from file name
% The regex captures patterns like 'energy_beam'_'beam_number'
% where 'energy' is the energy level, and 'beam_number' is the beam number
token = str2double(split(regexp(file_name, '\d+\.*\d*\_\d+', 'match', 'once'), "_"));
energy_beam = token(1);
beam_number = token(2);

%% Imports Data and check the file
% Opens the file to be processed
fide = fopen(file_name, 'r');
% Sorts each line in the .txt file into a cell array
data = textscan(fide, '%s', 'delimiter', '\n');
% Closes the file
fclose(fide);

% Finds the number of rows in the data file and stores it in n
n = size(data{1,1}, 1);
% Divides by 10 so that n represents the number of beam_numbers in the .txt file
n = n / (numDetect + 1);

% Checks that n is equal to the beam number. If not, an error is thrown.
if n ~= beam_number
    error('ERROR IN THE DATA FILE!');
end

%% Resets Counters
HitsLog = 0;
NoEnergyDep = 0;
percentage = 0;
index_count = 1;
index_arr = [];

%% Starts Post Process Loop
cd(outputfolder)
NewFileName = append('PostProcess', file_name, '.txt');
fid = fopen(NewFileName, 'wt');

fprintf('Starting %.2f \nPercent Complete:', energy_beam)
for i = 1:beam_number
    if rem(i, (beam_number / 10)) == 0
        percentage = percentage + 10;
        fprintf('%.0f ', percentage)
    end
    
    % Resets values for each iteration
    EnergySum = 0;
    reject_count = 0;
    
    % Preallocates Variables
    Detector_Energy = zeros(numDetect, 1);
    
    % Stores data from each line into Detector_Energy vector
    % Loop through all detectors: 1: " Edep (MeV)): "; 2-10: 9 detector readings
    % Line 10 ->Back Detector
    for j = 1:numDetect
        Detector_Energy(j) = str2double(data{1,1}{j + 1 + (i - 1) * (numDetect + 1), 1});
    end
    
    EnergySum = sum(Detector_Energy);
    
    if EnergySum == 0 % If no energy was deposited, do not transfer data to the output file
        NoEnergyDep = NoEnergyDep + 1;
    else % If energy is deposited, output hit to the output file
        % Writes Output to Text File
        fprintf(fid, ' Edep (MeV)): \n');
        fprintf(fid, '%.7f \n', Detector_Energy);
        HitsLog = HitsLog + 1;
        index_arr(index_count) = i; 
        index_count = index_count + 1;
    end
end

%% Writes Output to Text File
fprintf(fid, ' Sims with No Energy Deposited: \n');
fprintf(fid, '%.f \n', NoEnergyDep);
fprintf(fid, ' Sims with Energy Deposited: \n');
fprintf(fid, '%.f', (HitsLog));
fclose(fid);

% Return to MATLAB Main
cd(inputfolder);
cd ..
fprintf('\n %.2f Complete ', energy_beam);
end