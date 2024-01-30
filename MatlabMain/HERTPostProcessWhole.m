%% %% Authored by Skyler Krantz
% Edited by Will Teague 
% Updated: Nov. 13th, 2023
% Post Processing for GEANT4

function HERTPostProcessWhole(file_name, inputfolder, outputfolder)
% HERTPostProcess Removes simulations with no energy deposited
cd(inputfolder);

% Find the first sequence of digits and convert it to a number
beam_number = str2double(regexp(file_name, '\d+', 'match', 'once'));
run_number = str2double(regexp(file_name, '\d+_Run(\d+)', 'tokens', 'once'));

%% Imports Data and check the file
% Opens the file to be processed
fide = fopen(file_name, 'r');
% Sorts each line in the .txt file into a cell array
header = fgetl(fide);  % Read 11 strings directly into a cell array
data = textscan(fide, '%f %f %f %f %f %f %f %f %f %f', 'Delimiter','');  % Skip header
% Closes the file
fclose(fide);

% Separates data into incident energy and detected energy
Einc = data{1};
Detector_Energy = cell2mat(data(2:end));

% Finds the number of rows in the data file and stores it in n
n = size(data{1,1}, 1);

% Checks that n is equal to the beam number. If not, an error is thrown.
if n ~= beam_number
    error('ERROR IN THE DATA FILE!');
end

%% Starts Post Process Loop
cd(outputfolder)
NewFileName = append('PostProcess', file_name, '.txt');
fid = fopen(NewFileName, 'wt');

fprintf('Starting Run %d \n', run_number)
    
EnergySum = sum(Detector_Energy,2);

Einc_new = [];
Edep_data = [];

for i = 1:beam_number
    if EnergySum(i) ~= 0
        Einc_new = [Einc_new; Einc(i)];
        Edep_data = [Edep_data; Detector_Energy(i,:)];
    end
end

Energy_output = [Einc_new, Edep_data];

HitsLog = nnz(EnergySum);
NoEnergyDep = nnz(EnergySum == 0); 

%% Writes Output to Text File
fprintf(fid, '%s \n', header);
fprintf(fid, '%9.6g          %10.6g%10.6g%10.6g%10.6g%10.6g%10.6g%10.6g%10.6g%10.6g\n', Energy_output');
fprintf(fid, ' Sims with No Energy Deposited: \n');
fprintf(fid, '%.f \n', NoEnergyDep);
fprintf(fid, ' Sims with Energy Deposited: \n');
fprintf(fid, '%.f', HitsLog);
fclose(fid);

% Return to MATLAB Main
cd(inputfolder);
cd ..
fprintf('Run %d Complete \n', run_number);
end