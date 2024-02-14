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

n = size(Detector_Energy,1); % Number of particles in file

% Checks that n is equal to the beam number. If not, an error is thrown.
if n ~= beam_number
    error('ERROR IN THE DATA FILE!');
end

%% Starts Post Process Loop
cd(outputfolder)
NewFileName = append('PostProcess', file_name);
fid = fopen(NewFileName, 'wt');

fprintf('Starting Run %d \n', run_number)
    
EnergySum = sum(Detector_Energy,2);

% Use logical indexing and single assignments
idx_positive = EnergySum > 0;

Einc_new = Einc(idx_positive==1);
Einc_non = Einc(idx_positive==0);
Edep_data = Detector_Energy(idx_positive==1, :);

Energy_output = [Einc_new, Edep_data];

HitsLog = nnz(idx_positive);
NoEnergyDep = nnz(idx_positive == 0); 

%% Writes Output to Text File
fprintf(fid, 'Sims with Energy Deposited: %.f\n', HitsLog);
fprintf(fid, '%s \n', header);
fprintf(fid, '%9.6g          %10.6g%10.6g%10.6g%10.6g%10.6g%10.6g%10.6g%10.6g%10.6g\n', Energy_output');
fprintf(fid, '\nSims with No Energy Deposited: %.f\n',NoEnergyDep);
fprintf(fid, 'Einc(MeV) \n');
fprintf(fid, '%9.6g \n', Einc_non);
fclose(fid);

% Return to MATLAB Main
cd(inputfolder);
cd ..
fprintf('Run %d Complete \n', run_number);
end