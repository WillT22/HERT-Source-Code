function [hit_deposited_energy, hit_energy_channels, hit_detectors, energy_beam, run_number, beam_number, ...
    back_deposited_energy, back_energy_channels, back_detectors, back_energy_beam, non_energy_beam]...
    = oneEnergyEffDistWhole(filename, inputfolder, energy_channels, nucleons, detector_threshold, back_threshold)
% Author: Yinbo Chen
% Date: 6/15/2021
% Modified by: Skyler Krantz, Will Teague
% Date: Nov. 13, 2023

% Find the first sequence of digits and convert it to a number
beam_number = str2double(regexp(filename, '\d+', 'match', 'once'));
run_number = str2double(regexp(filename, 'Run(\d+)', 'tokens', 'once'));

cd(inputfolder);
fide = fopen(filename, 'r');
NumEnergyDeposit = textscan(fide, 'Sims with Energy Deposited: %f', 'Delimiter','');
deposit_data = textscan(fide, '%f %f %f %f %f %f %f %f %f %f', 'Delimiter','','HeaderLines',1); 
NumNoEnergy = textscan(fide, 'Sims with No Energy Deposited: %f', 'Delimiter','');
Einc_data = textscan(fide, '%f', 'Delimiter','','HeaderLines',1);  
fclose(fide);

NumEnergyDeposit = NumEnergyDeposit{1, 1};
energy_beam = deposit_data{1}'./nucleons;
Detector_Energy = cell2mat(deposit_data(2:end))./nucleons;
NumNoEnergy = NumNoEnergy{1, 1};
non_energy_beam = Einc_data{1}'./nucleons;

if NumNoEnergy + NumEnergyDeposit ~= beam_number
    error('ERROR IN THE DATA FILE!');
end
fprintf('Starting Run %d \n', run_number)

% 1. Evaluate Logic Vectors BEFORE Zeroing Data
count_reject_logic = Detector_Energy(:,1) < detector_threshold;         
back_hits_logic = Detector_Energy(:,end) > back_threshold;              

% 2. Apply Global Threshold
Detector_Energy(Detector_Energy < detector_threshold) = 0;

% 3. Process Rejected Hits
count_reject_indices = find(count_reject_logic);                        
non_energy_beam = [non_energy_beam, energy_beam(count_reject_indices)];  
count_reject = nnz(count_reject_logic);

% 4. Process Back Hits
% Define a single, strict mask for valid back hits
valid_back_hits_logic = back_hits_logic & ~count_reject_logic;
% Apply the exact same mask to all calculations
back_hits_indices = find(valid_back_hits_logic);                              
back_hits = nnz(valid_back_hits_logic);
back_energy_beam = energy_beam(valid_back_hits_logic);   
back_deposited_energy = sum(Detector_Energy(valid_back_hits_logic,:), 2);     
back_detectors = sum(Detector_Energy(valid_back_hits_logic,:) > 0, 1);

% 5. Process Good Hits (Exclude both rejected and back hits)
good_hits_logic = ~count_reject_logic & ~back_hits_logic;
energy_beam = energy_beam(good_hits_logic);                            
deposited_energy = Detector_Energy(good_hits_logic,:);                  
hit_detectors = sum(deposited_energy > 0, 1);             
hit_deposited_energy = sum(deposited_energy, 2); 

% 6. Assign Energy Channels for Good Hits
hit_energy_channels = zeros(1, length(energy_beam));
for ec = 1:size(energy_channels,1)                                                     
    for i = 1:size(deposited_energy,1)                                                       
        if hit_deposited_energy(i) >= energy_channels(ec, 1) && hit_deposited_energy(i) < energy_channels(ec, 2)      
            hit_energy_channels(i) = ec; 
        end
    end
end

% 7. Assign Energy Channels for Back Hits
back_energy_channels = zeros(1, length(back_energy_beam));
for ec = 1:size(energy_channels,1)                                                     
    for i = 1:size(back_deposited_energy,1)                                                       
        if back_deposited_energy(i) >= energy_channels(ec, 1) && back_deposited_energy(i) < energy_channels(ec, 2)      
            back_energy_channels(i) = ec; 
        end
    end
end

fprintf('Number of back hits= %i\n', back_hits);
fprintf('Number of rejected hits = %i\n', count_reject);
fprintf('Number of counted hits = %i\n', length(hit_energy_channels));

cd ..
end