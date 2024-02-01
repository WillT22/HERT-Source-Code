function [singleMatrix_whole, energy_beam, beam_number, count_back_whole, hits_detectors_whole, count_reject,run_number] = oneEnergyEffDistWhole(file_name, energy_channels, back_limit, inputfolder, detector_threshold)
% Author: Yinbo Chen
% Date: 6/15/2021
% Modified by: Skyler Krantz, Will Teague
% Date: Nov. 13, 2023

numDetect = 9;

% Initialize output variables
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

% Extract relevant information from the data
energy_beam = data{1};
Detector_Energy = cell2mat(data(2:end));
NumEnergyDeposit = NumEnergyDeposit{1, 1};
NumNoEnergy = NumNoEnergy{1, 1};

% Check if the sum of energy deposits and no-energy counts match the beam number
if NumNoEnergy + NumEnergyDeposit ~= beam_number
    error('ERROR IN THE DATA FILE!');
end

fprintf('Starting Run %d \n', run_number)

% Reset the counts
count_back_whole = 0;
count_reject = 0;

% Zero out values below detector threshold for the whole configuration
for i = 1:size(Detector_Energy,1)                       % For each particle that deposited nonzero energy on at least one detector,
    if Detector_Energy(i,1) < detector_threshold        % If the particle did not deposit above the threshold on the first detector,
        Detector_Energy(i,:) = 0;                       % Reject the count.
        count_reject = count_reject + 1;    
    elseif Detector_Energy(i,numDetect) > back_limit    % If the particle deposited energy above the threshold on the back detector,
        Detector_Energy(i,:) = 0;                       % Reject the count.
        count_back_whole = count_back_whole + 1;
    else                                                % Otherwise,
        for j = 2:numDetect-1                           % For the other detectors,
            if Detector_Energy(j) < detector_threshold  % If the particle did not deposit energy above the threshold of a detector,
                Detector_Energy(j) = 0;                 % set that detector's energy count to zero for that particle
            end
        end
    end
end
count_reject = count_reject + count_back_whole;

% Update hits count for each detector
for k = 1:(numDetect)
    hits_detectors_whole(k) = nnz(Detector_Energy(:,k));
end
    
% Calculate the sum of energy deposits over all detectors for each particle
WholeSum = sum(Detector_Energy,2); 

singleMatrix_whole = zeros(length(energy_beam), 1);

% Update singleMatrix for each energy channel
for l = 1:length(energy_channels(:, 1))                                                     % For each energy channel
    for i = 1:size(Detector_Energy,1)                                                       % For each particle
        if WholeSum(i) >= energy_channels(l, 1) && WholeSum(i) < energy_channels(l, 2)      % Sort total deposited energies into an energy channel
            singleMatrix_whole(i) = l;
        end
    end
end

% Display summary information after running through all simulations
fprintf('Number of back hits= %i\nTotal number of hits = %i\n', count_back_whole,sum(singleMatrix_whole));

% Change back to the original directory
cd ..

end