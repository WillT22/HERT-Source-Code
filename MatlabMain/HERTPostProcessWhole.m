function [Output] = HERTPostProcessWhole(file_name,inputfolder,outputfolder)
%HERTPostProcess Removes simulations with no energy deposited
numDetect = 9;

cd(inputfolder);

%recognize the paticle energy level from file name, only works when file
%name has 'energy_beam'_'beam_number' format number.
%\d+\. 1 or more digits followed by a period
%*\d*\_ any charcters plus 0 or more digits followed by a _. 
%Must use zero(*) in the case of whole number energy values. 
%This will capture energy_particlecount as well as energy formats
%\d+ 1 or more remaining digits
% Ex. Energy of 2 and beam number of 30 would result in 2_30
token = str2double(split(regexp(file_name, '\d+\.*\d*\_\d+', 'match', 'once'),"_"));
energy_beam = token(1);
beam_number = token(2);


%% Imports Data and check the file
n=0;
%Opens next file to be run
fide = fopen(file_name,'r');
%Sorts each line in the .txt file into a cell array
data = textscan(fide,'%s','delimiter','\n');

%Cloes files
fclose(fide);

%Finds in the number of rows in the data file and stores it to n
n = size(data{1,1},1);
%Divides by 10 so that n represents the number of beam_numbers in the .txt
%file
n = n/(numDetect+1);
%Checks that n is equal to the beam number. If not will stop the porgram
%with an error showing that something is wrong with the data file
if n ~= beam_number
    error('ERROR IN THE DATA FILE!');
end

%% Resets Counters
HitsLog=1;
NoEnergyDep = 0;
percentage = 0;

%% Starts Post Process Loop
cd(outputfolder)
NewFileName = append('PostProcess',file_name,'.txt');
fid = fopen(NewFileName,'wt');

fprintf('Starting %.2f \nPercent Complete:',energy_beam)
for i = 1:beam_number
   
    if rem(i,(beam_number/10)) == 0
        percentage = percentage + 10;
        fprintf('%.0f ',percentage)
    end
    
    %Resets values for each iteration
    EnergySum =0;
    reject_count=0;
    
    %Preallocates Variables
    Detector_Energy = zeros(numDetect,1);
    
    %Stores data from each line into Detector_Energy vector
    %Loop through all detectors: 1: " Edep (MeV)): "; 2-10: 9 detector readings
    %Line 10 ->Back Detector
    for j = 1:numDetect
        Detector_Energy(j) = str2double(data{1,1}{j+1+(i-1)*(numDetect+1),1});
        
    end
    
    EnergySum = sum(Detector_Energy);
    %if no energy was deposited, do not transfer data to output file
    if EnergySum == 0
        NoEnergyDep = NoEnergyDep +1;
    else
    %if energy is deposited, output hit to output file
        % Writes Output to Text File
        fprintf(fid,' Edep (MeV)): \n');
        fprintf(fid,'%.7f \n',Detector_Energy);
        HitsLog= HitsLog+1;
    end
end

%% Writes Output to Text File
fprintf(fid,' Sims with No Energy Deposited: \n');
fprintf(fid,'%.f \n',NoEnergyDep );
fprintf(fid,' Sims with Energy Deposited: \n');
fprintf(fid,'%.f',(HitsLog-1));
fclose(fid);

%Lets loop in Main0_HERTPostProcessLoop.m move to next count
Output =1;

%Return to MATLAB Main

cd(inputfolder);
cd ..
fprintf('\n %.2f Complete \n',energy_beam);

end

