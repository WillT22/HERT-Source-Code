function [singleMatrix_whole,energy_beam,beam_number,hits_log,count_back_whole,detector_energy_whole,hits_detectors_whole] = oneEnergyEffDistWhole(file_name,energy_channels,outer_limit,inputfolder,detector_threshold)
%Author:Yinbo Chen
%Date: 6/15/2021
%Modified by: Skyler Krantz
%Date 6/15/2022
numDetect = 9;
% %Begins Timer
% tic
%input txt
cd(inputfolder);

%initial output x*1 array with name singleMatrix

singleMatrix_whole = zeros(length(energy_channels), 1);
hits_log = 0;


hits_detectors_whole= zeros(numDetect,1);

detector_energy_whole= zeros(numDetect,1);

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
%number = zeros(length(beam_number),18);

%% Imports Data and check the file
n=0;
%Opens next file to be run
fide = fopen(file_name,'r');
%Sorts each line in the .txt file into a cell array
data = textscan(fide,'%s','delimiter','\n');

%Cloes files
fclose(fide);


%Finds in the number of rows in the datat file and stores it to n
[n,~] = size(data{1,1});
%Divides by 18 so that n represents the number of beam_numbers in the .txt
%file

NumEnergyDeposit = str2double(data{1,1}{n,1});
NumNoEnergy = str2double(data{1,1}{n-2,1});
%Checks that n is equal to the beam number. If not will stop the porgram
%with an error showing that something is wrong with the data file
if NumNoEnergy+NumEnergyDeposit ~= beam_number
    error('ERROR IN THE DATA FILE!');
end


%% Resets Counter
count_outer_combined = 0;
count_back_whole = 0;


percentage = 0;
fprintf('Starting %.2f \nPercent Complete:',energy_beam)
h = 1;
s = 1;

%% Loop through all simulations in filename
for i = 1:NumEnergyDeposit
    %Prints Percentage complete of the file
    if i>(percentage*(10^-2)*NumEnergyDeposit)
        percentage = percentage + 10;
        fprintf('%.0f ',percentage)
    end
    
    %Resets values for each iteration
    
    
    WholeSum = 0;
    reject_count=0;
    
    %Preallocates Variables
    Detector_Energy = zeros(numDetect,1);
    Detector_Energy = zeros(numDetect,1);
    
    %Stores data from each line into Detector_Energy vector
    %Loop through all detectors: 1: " Edep (MeV)): "; 2-18: 17 detector readings
    %Lines 2,4,6,8,10,12,14,16 -> Inner Detectors
    %Lines 3,5,7,9,11,13,15,17 -> Outer Detectors
    %Line 18 ->Back Detector
    for j = 1:numDetect
        Detector_Energy(j) = str2double(data{1,1}{j+1+(i-1)*(numDetect+1),1});
        if Detector_Energy(j)>0
            hits_detectors_whole(j) = hits_detectors_whole(j)+1;
        end
    end
    
    %     Whole_Detector_Energy(1) = Detector_Energy(1);
    %     %Whole_Detector_Energy(2) = Detector_Energy(3);
    %     if Whole_Detector_Energy(1)>0
    %         hits_detectors_whole(1) = hits_detectors_whole(1)+1;
    %     end
        
%     for k = 2:9
%         if k <= 8
%         Whole_Detector_Energy(k) = Detector_Energy((k-1)*2+1) + Detector_Energy((k-1)*2+2);
%         else
%             Whole_Detector_Energy(k)=Detector_Energy(17);
%         end
%         % If any energy is deposited, records a hit on that detector
%         if Whole_Detector_Energy(k)>0
%             hits_detectors_whole(k) = hits_detectors_whole(k)+1;
%         end
%         
% 
%     end
    
    
    
    % detector_energy_inner = detector_energy_inner + Detector_Energy;
    detector_energy_whole = detector_energy_whole + [Detector_Energy];
    
    %Checks if the back detector has energy greater than threshold.
    %If so, it skips to next simulation
    %Changed to 100 for side pen tests
    if Detector_Energy(numDetect) > outer_limit
        count_back_whole = count_back_whole +1;
        WholeSum = 0;
        
        reject_count=1;
        continue
    end
    
    %For whole configuration, zeros out values below detector threshold 
    for k = 1:numDetect
        if Detector_Energy(k) < detector_threshold
            Detector_Energy(k) = 0;
        end
    end
    
%     %For inner configuration, zeros out detectors below threshold
%     for L = 1:16
%         if Detector_Energy(L) < detector_threshold
%             Detector_Energy(L) = 0;
%         end
%     end

    %Creates Sums for Both Configurations
    WholeSum = sum(Detector_Energy);
%     InnerSum = sum(Detector_Energy(1:2:15));
%     OuterSum = sum(Detector_Energy(2:2:16));
    
%     %Filters out Hit for inner ring configuration if it exceeds outer_limit
%     if OuterSum > outer_limit
%         count_outer_inner = count_outer_inner+1;
%         InnerSum = 0;
%     end
       
    %Loop through all energy channels. If the energy deposited in the
    %detectors is in between the energy channels it will add a count to
    %that energy channel

    if WholeSum > 0  && Detector_Energy(1) > detector_threshold && reject_count == 0
        for k = 1:length(energy_channels(:,1))
            if WholeSum >= energy_channels(k,1) && WholeSum < energy_channels(k,2)
                singleMatrix_whole(k) = singleMatrix_whole(k)+1;
                hits_log(h) = i;
                h= h+1;
                
                %detector_energy_whole = detector_energy_whole + [Whole_Detector_Energy(:);Detector_Energy(17)];
            end
        end
    end
    
    %     if InnerSum > 0 && Detector_Energy(1)  > detector_threshold && reject_count == 0
    %         for k = 1:length(energy_channels(:,1))
    %             if InnerSum >= energy_channels(k,1) && InnerSum < energy_channels(k,2)
    %                 singleMatrix_inner(k) = singleMatrix_inner(k)+1;
    %                 hits_log_inner(s) = i;
    %                 s= s+1;
    %
    %                 %detector_energy_inner = detector_energy_inner + Detector_Energy;
    %             end
    %         end
    %     end
    
end
%After running through all simulations in the file, singleMatrix will
%contain how many simulations showed energy in each energy channel.
%Ex. for 10000 simulations, if 3000 were in energy channel 1
%without energy in the outer_ring or hitting the back
%then singleMatrix(1) = 3000

%Divides deposited energy by the number of hits to get average energy
%deposited per particle simulated
%detector_energy_inner = detector_energy_inner./beam_number;
%detector_energy_whole = detector_energy_whole./beam_number;

%Prints the current energy level and how many particles hit the back wall
%or detector
fprintf('\n Energy Level:%.2i\n',energy_beam)
fprintf('Whole Configuration: The number of beams hit the outer_ring=%i\nThe numbers of beams hit the back=%i\nThe numbers of beams hit =%i\n',count_outer_combined,count_back_whole,sum(singleMatrix_whole));
% fprintf('Inner Ring: The number of beams hit the outer_ring=%i\nThe numbers of beams hit the back=%i\nThe numbers of beams hit =%i\n',count_outer_inner,count_back_inner,sum(singleMatrix_inner));


%Divides singleMatrix by how many simulations in the file to convert to a
%percentage. Ex. for 10000 simulations if 3000 were in energy channel 1,
%then is now singleMatrix(1) = 3000/10000 = 0.3
%singleMatrix = singleMatrix/beam_number;

cd ..

% %Ends Timer
% toc
end