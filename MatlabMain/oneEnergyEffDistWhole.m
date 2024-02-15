function [hit_energy_channels, run_number, beam_number, energy_beam, non_energy_beam, back_energy_beam] = oneEnergyEffDistWhole(file_name, inputfolder, energy_channels, detector_threshold, back_threshold)
% Author: Yinbo Chen
% Date: 6/15/2021
% Modified by: Skyler Krantz, Will Teague
% Date: Nov. 13, 2023

% Find the first sequence of digits and convert it to a number
beam_number = str2double(regexp(file_name, '\d+', 'match', 'once'));
run_number = str2double(regexp(file_name, '\d+_Run(\d+)', 'tokens', 'once'));

% Change to the input folder
cd(inputfolder);

% Opens the file to be processed
fide = fopen(file_name, 'r');
% Sorts each line in the .txt file into a cell array
NumEnergyDeposit = textscan(fide, 'Sims with Energy Deposited: %f', 'Delimiter','');
deposit_data = textscan(fide, '%f %f %f %f %f %f %f %f %f %f', 'Delimiter','','HeaderLines',1);  % Skip header
NumNoEnergy = textscan(fide, 'Sims with No Energy Deposited: %f', 'Delimiter','');
Einc_data = textscan(fide, '%f', 'Delimiter','','HeaderLines',1);  % Skip header
% Closes the file
fclose(fide);

% Extract relevant information from the data
NumEnergyDeposit = NumEnergyDeposit{1, 1};
energy_beam = deposit_data{1}';
Detector_Energy = cell2mat(deposit_data(2:end));
NumNoEnergy = NumNoEnergy{1, 1};
non_energy_beam = Einc_data{1}';

% Check if the sum of energy deposits and no-energy counts match the beam number
if NumNoEnergy + NumEnergyDeposit ~= beam_number
    error('ERROR IN THE DATA FILE!');
end

fprintf('Starting Run %d \n', run_number)

% Zero out values below detector threshold for the whole configuration
Detector_Energy(Detector_Energy < detector_threshold) = 0;

% Counts rejected hits
count_reject_logic = Detector_Energy(:,1)==0;                           % find which hits are rejected
count_reject_indices = find(count_reject_logic);                        % find indices of rejected hits
non_energy_beam = [non_energy_beam,energy_beam(count_reject_indices)];  % update non_energy_beam with rejected hits
energy_beam = energy_beam(~count_reject_logic);                         % remove energy values from energy beam
Detector_Energy = Detector_Energy(~count_reject_logic,:);               % update Detector_Energy with non-rejected hits

% Counts back hits
back_energy_beam = [];
back_hits_logic = Detector_Energy(:,end)> back_threshold;               % find which hits are rejected
back_hits_indices = find(back_hits_logic);                             % find indices of rejected hits
back_energy_beam = [back_energy_beam,energy_beam(back_hits_indices)];   % update non_energy_beam with rejected hits
energy_beam = energy_beam(~back_hits_logic);                             % remove energy values from energy beam
Detector_Energy = Detector_Energy(~back_hits_logic,:);                  % update Detector_Energy with non-rejected hits

back_hits = nnz(back_hits_logic);
count_reject = nnz(count_reject_logic) + nnz(back_hits_logic);
    
% Calculate the sum of energy deposits over all detectors for each particle
WholeSum = sum(Detector_Energy,2); 

hit_energy_channels = zeros(1,length(energy_beam));
% Update singleMatrix for each energy channel
for ec = 1:size(energy_channels,1)                                                     % For each energy channel
    for i = 1:size(Detector_Energy,1)                                                       % For each particle
        if WholeSum(i) >= energy_channels(ec, 1) && WholeSum(i) < energy_channels(ec, 2)      % Sort total deposited energies into an energy channel
            hit_energy_channels(i) = ec; %replace with variable sized matrix that takes all good hits
        end
    end
end

% Display summary information after running through all simulations
fprintf('Number of back hits= %i\n', back_hits);
fprintf('Number of rejected hits = %i\n', count_reject);
fprintf('Number of counted hits = %i\n', length(hit_energy_channels));

% Change back to the original directory
cd ..

end