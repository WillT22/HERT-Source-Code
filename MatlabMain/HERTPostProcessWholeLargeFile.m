function [Output] = HERTPostProcessWholeLargeFile(file_name,inputfolder,outputfolder)
%HERTPostProcess Removes simulations with no energy deposited
%Large File reads the file line by line. This takes more time but limits
%memory usage.
% ------------------THIS CODE IS NOT COMPUTATIONALLY OPTIMIZED------------

numDetect = 9;

cd(inputfolder);

%recognize the paticle energy level from file name, only works when file
%name has 3.8_10 format number.
%3.8 is energy; 10 is beam number
%\d+\. 1 or more digits followed by a period
%*\d*\_ any charcters plus 0 or more digits followed by a _. Must use zero(*) in the case of
%whole number energy values. This will capture 3.8_10 as well as 3.30
%\d+ 1 or more remaining digits
% Ex. Energy of 2 and beam number of 30 would result in 2_30
token = str2double(split(regexp(file_name, '\d+\.*\d*\_\d+', 'match', 'once'),"_"));
energy_beam = token(1);
beam_number = token(2);


%% Imports Data and check the file
n=0;
%Opens next file to be run
infile = fopen(file_name,'r');

%% Commenting out for large files
% %Sorts each line in the .txt file into a cell array
% data = textscan(infile,'%s','delimiter','\n');
% 
% %Cloes files
% %fclose(infile);
% 
% %Finds in the number of rows in the datat file and stores it to n
% [n,~] = size(data{1,1});
% %Divides by 18 so that n represents the number of beam_numbers in the .txt
% %file
% n = n/(numDetect+1);
% %Checks that n is equal to the beam number. If not will stop the porgram
% %with an error showing that something is wrong with the data file
% if n ~= beam_number
%     error('ERROR IN THE DATA FILE!');
% end

%% Resets Counters
HitsLog=1;
NoEnergyDep = 0;
percentage = 0;

%% Starts Post Process Loop
cd(outputfolder)
NewFileName = append('PostProcess',file_name,'.txt');
outfile = fopen(NewFileName,'wt');

fprintf('Starting %.2f \nPercent Complete:',energy_beam)
for i = 1:beam_number
    Detector_Energy = zeros(numDetect,1);
    % Reads in next line and stores to data
    % Stores detector energy for each simulation
    for line = 1:(numDetect+1)
        if line>1
        Detector_Energy(line-1) = str2double(fgetl(infile));%reads energy deposited line
        else 
            fgetl(infile);%reads header line and goes to next line
        end
        
    end
    
    if rem(i,(beam_number/10)) == 0
        percentage = percentage + 10;
        fprintf('%.0f ',percentage)
    end
    
    
    %Resets values for each iteration
    EnergySum =0;
    
    reject_count=0;
        
    %Stores data from each line into Detector_Energy vector
    %Loop through all detectors: 1: " Edep (MeV)): "; 2-18: 17 detector readings
    %Lines 2,4,6,8,10,12,14,16 -> Inner Detectors
    %Lines 3,5,7,9,11,13,15,17 -> Outer Detectors
    %Line 18 ->Back Detector
%     for j = 1:numDetect
%         Detector_Energy(j) = str2double(data{1,1}{j+1+(i-1)*(numDetect+1),1});
%         
%     end
    
    EnergySum = sum(Detector_Energy);
    if EnergySum == 0
        NoEnergyDep = NoEnergyDep +1;
    else
        
        % Writes Output to Text File
        fprintf(outfile,'Edep (MeV)): \n');
        fprintf(outfile,'%.7f \n',Detector_Energy);
        HitsLog= HitsLog+1;
    end
    
    
end

fclose(infile);
%% Writes Output to Text File
fprintf(outfile,' Sims with No Energy Deposited: \n');
fprintf(outfile,'%.f \n',NoEnergyDep );
fprintf(outfile,' Sims with Energy Deposited: \n');
fprintf(outfile,'%.f',(HitsLog-1));
fclose(outfile);

Output =1;

%Return to MATLAB Main\Result
cd 'D:\HERT\MATLAB Main\Result'
fprintf('\n %.2f Complete \n',energy_beam);

end

