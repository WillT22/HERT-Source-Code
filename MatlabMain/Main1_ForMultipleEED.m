%%HERT: mainForMultipleEED.mfiletype
% Last modified: 5/8/2023


% This script is the main code for converting Geant4 results to derived
% count rate and geometric factor
% Requires MatPlotLib Perceptually Uniform Colormaps

%Resets all variables and values in MATLAB
clear;
close all;
clc;
addpath 'E:\HERT_Drive\MATLAB Main'

%% Geometric Factor- Theory
% This section solves for the theorectical geometric factor of the instrument
% assuming the collimator teeth are perfect.
% Calculation dervied from Sullivan 1971 paper:
% https://www.sciencedirect.com/science/article/abs/pii/0029554X71900334

%% Assuming the Collimator Knife Edge Stops All Particles (Min GeoFactor)
% Disc 1: First Collimator Tooth
% Disc 2: Last Collimator Tooth
% Disc 3: First Detector

L_12 = 6.0;%cm distance between first and last collimator teeth
L_23 = 0.3; %cm distance between last collimator tooth and first detector
r1 = 0.9; %cm radius of first collimator tooth
r2 = 0.9; %cm radius of last collimator tooth
r3 = 1.0; %cm radius of first detector

L_13 = L_12+L_23;

%Creating angles for simplifying criteria Eq. 14 in Sullivan
%Angles of incidence such that circles project onto each other depending on location of incidence of the particles
%maximum angle that the particle will always hit within both circles
theta_c_12 = atan(abs(r1-r2)/L_12); %first to last collimator tooth
theta_c_13 = atan(abs(r1-r3)/L_13); %first tooth to first detector
theta_c_23 = atan(abs(r2-r3)/L_23); %last tooth to detector

%maximum angle such that particles may hit both circles depending on incidence location
theta_m_12 = atan((r1+r2)/L_12); %first to last collimator tooth
theta_m_13 = atan((r1+r3)/L_13); %first tooth to first detector
theta_m_23 = atan((r2+r3)/L_23); %last tooth to detector

%geometric factor of first collimator tooth and first detector (high E particles)
G13 = 0.5*(pi^2)*((r1^2+r3^2+L_13^2)-(((r1^2+r3^2+L_13^2)^2-4*(r1^2)*(r3^2))^0.5));
%geometric factor of last collimator tooth and first detector (low E particles)
G12 = 0.5*(pi^2)*((r1^2+r2^2+L_12^2)-((r1^2+r2^2+L_12^2)^2-4*(r1^2)*(r2^2))^0.5); 
%geometric factor of first to last collimator tooth
G23 = 0.5*(pi^2)*((r2^2+r3^2+L_23^2)-((r2^2+r3^2+L_23^2)^2-4*(r2^2)*(r3^2))^0.5); 

%Apply Simlifying Criteria
%if the 'always hits' critical angle is defined from the first tooth and the detector
if theta_c_12 >= theta_c_13
    G13 = 0.5*(pi^2)*((r1^2+r3^2+L_13^2)-(((r1^2+r3^2+L_13^2)^2-4*(r1^2)*(r3^2))^0.5));
    G3_whole = G13; %geometric factor is defined by the first tooth and first detector
    fprintf('Geometric_Factor_13 (First Tooth to First Detector)= %7.5f  cm^2 sr\n \n',G13)
    FOV = 2*theta_m_13*180/pi;

%otherwise, if the 'can hit' critical angle is defined by the collimator
elseif theta_m_13 >= theta_m_12
    G12 = 0.5*(pi^2)*((r1^2+r2^2+L_12^2)-((r1^2+r2^2+L_12^2)^2-4*(r1^2)*(r2^2))^0.5);
    G3_whole = G12; %geometric factor is defined by the collimator
    fprintf('Geometric_Factor_12 (First Tooth to Last Tooth) = %7.5f  cm^2 sr\n \n',G12)
    FOV = 2*theta_m_12*180/pi;
    
%otherwise, if the 'can hit' critical angle is defined by the last tooth and the detector
elseif theta_c_12 >= theta_m_13
    G23 = 0.5*(pi^2)*((r2^2+r3^2+L_23^2)-((r2^2+r3^2+L_23^2)^2-4*(r2^2)*(r3^2))^0.5);
    G3_whole = G23; %geometric factor is defined by the last tooth and first detector
    fprintf('Geometric_Factor_23 (Last Tooth to First Detector) = %7.5f  cm^2 sr\n \n',G23)

%otherwise, the geometric factor is defined by all three components
else
    theta_a = atan(((L_23*r1^2+L_12*r3^2-L_13*r2^2)^0.5)/(L_12*L_23*L_13));
    G13 = 0.5*(pi^2)*((r1^2+r3^2+L_13^2)-(((r1^2+r3^2+L_13^2)^2-4*(r1^2)*(r3^2))^0.5));
    
    %from Sullivan eq 16
    Z12= Zij(theta_a,L_12,r1,r2,theta_c_12,theta_m_12,r2);
    Z13= Zij(theta_a,L_13,r1,r3,theta_c_13,theta_m_13,r2);
    Z23= Zij(theta_a,L_23,r2,r3,theta_c_23,theta_m_23,r2);
    
    G123=G13-pi^2*r2^2*sin(theta_a)^2+Z23+Z12-Z13; %Sullivan eq 15
    G3_whole = G123;
    fprintf('Geometric_Factor (Three Disc Telescope) = %7.5f  cm^2 sr\n \n',G123)
    
end

%% Assuming the Collimator Knife Edge Stops No Particles (Max GeoFactor)
% Disc 1: First Collimator Tooth
% Disc 2: Last Collimator Tooth
% Disc 3: First Detector

L_12 = 6.0;%cm distance between first and last collimator teeth
L_23 = 0.3; %cm distance between last collimator tooth and first detector
r1 = 1.0; %cm radius of first collimator tooth (larger than above)
r2 = 1.0; %cm radius of last collimator tooth (larger than above)
r3 = 1.0; %cm radius of first detector

L_13 = L_12+L_23;

theta_c_12 = atan(abs(r1-r2)/L_12);
theta_c_13 = atan(abs(r1-r3)/L_13);
theta_c_23 = atan(abs(r2-r3)/L_23);

theta_m_12 = atan((r1+r2)/L_12);
theta_m_13 = atan((r1+r3)/L_13);
theta_m_23 = atan((r2+r3)/L_23);

G13 = 0.5*(pi^2)*((r1^2+r3^2+L_13^2)-(((r1^2+r3^2+L_13^2)^2-4*(r1^2)*(r3^2))^0.5));
G12 = 0.5*(pi^2)*((r1^2+r2^2+L_12^2)-((r1^2+r2^2+L_12^2)^2-4*(r1^2)*(r2^2))^0.5);
G23 = 0.5*(pi^2)*((r2^2+r3^2+L_23^2)-((r2^2+r3^2+L_23^2)^2-4*(r2^2)*(r3^2))^0.5);

if theta_c_12 >= theta_c_13
    G13 = 0.5*(pi^2)*((r1^2+r3^2+L_13^2)-(((r1^2+r3^2+L_13^2)^2-4*(r1^2)*(r3^2))^0.5));
    G3_whole_max = G13;
    fprintf('Geometric_Factor_13 (Front Tooth to First Detector)= %7.5f  cm^2 sr\n \n',G13)
    FOV = 2*theta_m_13*180/pi
    
elseif theta_m_13 >= theta_m_12
    
    G12 = 0.5*(pi^2)*((r1^2+r2^2+L_12^2)-((r1^2+r2^2+L_12^2)^2-4*(r1^2)*(r2^2))^0.5);
    G3_whole_max = G12;
    fprintf('Geometric_Factor_12 (First Tooth to Last Tooth) = %7.5f  cm^2 sr\n \n',G12)
    FOV = 2*theta_m_12*180/pi
    
elseif theta_c_12 >= theta_m_13
    G23 = 0.5*(pi^2)*((r2^2+r3^2+L_23^2)-((r2^2+r3^2+L_23^2)^2-4*(r2^2)*(r3^2))^0.5);
    G3_whole_max = G23;
    fprintf('Geometric_Factor (Last Tooth to First Detector) = %7.5f  cm^2 sr\n \n',G23)
    
else
    theta_a = atan(((L_23*r1^2+L_12*r3^2-L_13*r2^2)^0.5)/(L_12*L_23*L_13));
    G13 = 0.5*(pi^2)*((r1^2+r3^2+L_13^2)-(((r1^2+r3^2+L_13^2)^2-4*(r1^2)*(r3^2))^0.5));
    
    Z12= Zij(theta_a,L_12,r1,r2,theta_c_12,theta_m_12,r2);
    Z13= Zij(theta_a,L_13,r1,r3,theta_c_13,theta_m_13,r2);
    Z23= Zij(theta_a,L_23,r2,r3,theta_c_23,theta_m_23,r2);
    
    G123=G13-pi^2*r2^2*sin(theta_a)^2+Z23+Z12-Z13;
    G3_whole_max = G123;
    fprintf('Geometric_Factor (Three Disc Telescope) = %7.5f  cm^2 sr\n \n',G123)
    
end

%% Inner/Outer
rI= zeros(1,8);
theta_c_inner= zeros(length(rI));
L_inner = zeros(length(rI));
Length_inner = 0.25; %cm
theta_m_inner= zeros(length(rI));
r_coll=1;%cm
l_coll=7.45;%cm
%
%Inner Ring Only
for i = 1:8
    rI(i)= 1; %cm
end
%
%
for i = 1:length(rI)
    for j = 1:length(rI)
        if i > j
            L_inner(i,j) = 0/0;
            theta_c_inner(i,j) = 0/0;
            theta_m_inner(i,j) = 0/0;
            
        else
            L_inner(i,j) = (j-i)*Length_inner;
            theta_m_inner(i,j) = atand((rI(i)+rI(j))/L_inner(i,j));
            theta_c_inner(i,j) = atand(abs((rI(i)-rI(j)))/L_inner(i,j));
        end
    end
end


for i = 1:length(rI)
    theta_c_coll(i) = atand(abs((r_coll-rI(i)))/(l_coll+L_inner(1,i)));
    theta_m_coll(i) = atand((r_coll+rI(i))/(l_coll+L_inner(1,i)));
end

if max(min(theta_c_inner,[],1)) == theta_c_inner(1,8)
    Gtest = 0.5*(pi^2)*((rI(8)^2+rI(1)^2+L_inner(1,8)^2)-(((rI(8)^2+rI(1)^2+L_inner(1,8)^2)^2-4*(rI(8)^2)*(rI(1)^2))^0.5));
    a_x = rI(8)^2+rI(1)^2+L_inner(1,8)^2;
    b_y = 4*(rI(8)^2)*(rI(1)^2);
    Gtest1_inner = 0.5*(pi^2)*(a_x-(a_x^2-b_y)^0.5);
end
%Gtest1_inner = Gtest *30/180
G1_inner = 0.5*(pi^2)*((r_coll^2+rI(1)^2+l_coll^2)-(((r_coll^2+rI(1)^2+l_coll^2)^2-4*(r_coll^2)*(rI(1)^2))^0.5))
G2_inner = 0.5*(pi^2)*((r_coll^2+rI(8)^2+(l_coll+L_inner(1,8))^2)-(((r_coll^2+rI(8)^2+(l_coll+L_inner(1,8))^2)^2-4*(r_coll^2)*(rI(8)^2))^0.5))
%G3_whole = 0.5*(pi^2)*((r_coll^2+rI(1)^2+(l_coll)^2)-((r_coll^2+rI(1)^2+(l_coll)^2)^2-4*(r_coll^2)*(rI(1)^2))^0.5);
%
%
% G3_test = 0.5*(pi^2)*((r_coll^2+r_coll^2+(l_coll)^2)-((r_coll^2+r_coll^2+(l_coll)^2)^2-4*(r_coll^2)*(r_coll^2))^0.5);
% G3_whole = G3_test;
% %fprintf('Geometric_Factor (Front Coll. to Detector 1) = %7.5f  cm^2 sr\n \n',G1_inner)
% %fprintf('Geometric_Factor_Inner (Front Coll. to Detector 15)= %7.5f  cm^2 sr\n \n',G2_inner)
% fprintf('Geometric_Factor_Whole (Front Coll. to First Detector)= %7.5f  cm^2 sr\n \n',G3_whole)
%
% fprintf('Geometric_Factor_Whole Test = %7.5f  cm^2 sr\n \n',G3_test)
%


%% Calculate GEANT4 Results
%read files from ./Result folder stores into 1*C array
cd 'E:\HERT_Drive\MATLAB Main\Result'; %Main Result Directory

% Get Folder Names for User
% (Source:https://www.mathworks.com/matlabcentral/answers/166629-is-there-...
% any-way-to-list-all-folders-only-in-the-level-directly-below-a-selected-directory)
topLevelFolder = pwd; % or whatever, such as 'C:\Users\John\Documents\MATLAB\work'
% Get a list of all files and folders in this folder.
files = dir(topLevelFolder);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags); % A structure with extra info.
% Get only the folder names into a cell array.
subFolderNames = {subFolders(3:end).name}; % Start at 3 to skip . and ..

FolderChoice = menu('HERT Loop: Choose an input folder',subFolderNames{:});
inputfolder = subFolderNames{FolderChoice};
addin = inputfolder;
%Menu to select spherical cap or full spherical
simtype_choice = menu('Simulation Type','Spherical Cap (15 deg)','Full Sphere');
switch simtype_choice
    
    case 1
        sim_type = 0;
        fprintf('Sim Type: Spherical Cap (15 deg) \n')
        addin = append(addin,' SC ');
        
    case 2
        sim_type = 1;
        fprintf('Sim Type: Full Sphere \n')
        addin = append(addin,' FS ');
        
end

%Creates search string for result .txt files
inputfiles = append(inputfolder,'/*.txt');

%lists all .txt files in the Result folder
list = dir(inputfiles);

%Grabs all the names of the files in vector (nx1 matrix)
list_fileNames = {list.name};

%Gets number of rows and columns in file names. Columns will indicate
%number of files that can be loaded
[R,C] = size(list_fileNames);

% Display a menu and get a choice
choice = menu('Choose an option', 'Exit Program', 'Load one file','Load all files');
%Exit Program =1 Load one file = 2 Load all files =3 Start Run=4

% Choice 1 is to exit the program
while choice ~= 1
    switch choice
        case 0
            disp('Error - please choose one of the options.')
            
            %Load One File
        case 2
            %Displays all the files that can be loaded
            txt_file_choice = menu('Choose a file',list_fileNames{:});
            
            %If menu button is closed, it will recycle to inital menu
            if txt_file_choice == 0
                disp('Please select a file')
                %If user selects a file, it will change filename to reference that one file.
            elseif txt_file_choice > 0
                for i = 1:C
                    switch txt_file_choice
                        case i
                            filename = list_fileNames{i};
                    end
                end
                %Shows the file that was loaded
                disp(filename);
            end
            
            %Load all files
        case 3
            txt_files_All= menu('Choose all files','All .txt files');
            switch txt_files_All
                case 1
                    %Loads all .txt files in Result folder
                    filename = list_fileNames';
                    fprintf('Number of files loaded: %.0f',length(list_fileNames))
                    
            end
            
        % Start Run
        case 4
            %Menu to determine runtype
            runtype_choice = menu('Outer Ring and Back Detector Limit','0 MeV','.1 MeV limit','0.01 MeV limit');
            switch runtype_choice
                %Includes Outer Ring Hits
                case 1
                    outer_limit = 0;
                    addin = append(' 0 MeV OTR ',addin);
                    %Excludes Outer Ring Hits
                case 2
                    outer_limit = 0.1;
                    addin = append(' 0.1 MeV OTR ',addin);
                    %Sets 0.1 MeV Threshhold for Outer Ring
                case 3
                    outer_limit = 0.01;
                    addin = append(' 0.01 MeV OTR ',addin);
            end
            %Prints Outer_limit
            fprintf('Outer_limit:%.2i\n',outer_limit)
            
            %Menu to determine detector threshold
            detector_choice = menu('Detector Threshold','0 MeV','0.1 MeV');
            switch detector_choice
                
                case 1
                    detector_threshold = 0;
                    addin = append(addin,' 0 MeV DT ');
                    
                case 2
                    detector_threshold = 0.1; %MeV
                    addin = append(addin,' 0.1 MeV DT ');
            end
            %Prints Outer_limit
            fprintf('Detector Threshold:%.2i MeV\n',detector_threshold)
            
            %Menu to determine energy channels
            %Finds all .txt files in Results directory. These txt files
            %contain the energy channels
            channels = dir('*.txt');
            
            %Produces menu for user to select while energy channels
            energy_channel_choice = menu('Choose Energy Channels',channels.name);
            
            %Reads in selected energy channel file and prepares addin for result
            %plots
            energy_channels = readmatrix(channels(energy_channel_choice).name);
            Selected_Channel_name = channels(energy_channel_choice).name;
            size_EC = size(energy_channels);
            fprintf('%f Energy Channels Selected \n',size_EC(1))
            addin = append(' ',erase(channels(energy_channel_choice).name,'.txt'),addin);
            
            %Preallocates EngLegend variable
            EngLegend = strings([1,length(energy_channels)]);
            
            %Creates Str Array for Plot legend
            for i = 1:size_EC(1)
                EngLegend(i) = append(num2str(round(energy_channels(i,1),2)),'-',num2str(round(energy_channels(i,2),2)),' MeV');
            end
            
            %Creates n x 3 matrix for plot colors. This will give each energy channel its own color on the plot.
            %n = number of energy channels
            Effplotcolor = plasma(size_EC(1)); %requires MatPlotLib Perceptually Uniform Colormaps
            
            %One file selected
            if size(filename,1)==1
                disp('Start the oneEnergyEffDist.m');  
                %Runs oneEnergyEffDistWhole for the one .txt file
                [output_Mult,output_energy,output_number,hits_log,count_back_whole,detector_energy_whole,hits_detectors_whole]...
                    = oneEnergyEffDistWhole(filename,energy_channels,outer_limit,inputfolder,detector_threshold);
                    
                hits_whole = sum(output0.*hits_whole);
                
                save('output_singleParticleArray.mat','output0')
                disp('output_singleParticleArray.mat');
                x= linspace(1.0,7.0,length(energy_channels));
                y_whole = output0' ;
                figure
                plot(x,y_whole(:,:))
                
                % More than one file selected
            elseif size(filename,1) > 1
                disp('Start to loop oneEnergyEffDist.m');
                disp(addin);
                
                %Creates matrix to store data
                final_Matrix = zeros(length(energy_channels),C);
                final_Matrix2 = zeros(length(energy_channels),C);
                final_Matrix3 = zeros(length(energy_channels),C);
                final_Matrix4 = zeros(length(energy_channels),C);
                final_Matrix5 = zeros(1,C);
                final_Matrix6 = zeros(1,C);
                final_Matrix7 = zeros(1,C);
                final_Matrix8 = zeros(9,C);
                final_Matrix9 = zeros(17,C);
                final_Matrix10 = zeros(9,C);
                final_Matrix11 = zeros(17,C);
                
                %Nested For loops to create final matrix 1 and 2
                for i = 1: C
                    tic
                    %For every .txt file in Results, it will run
                    %oneEnergyEffDist and add the results to finalMatrix
                    %and finalMatrix2
                    [output_Mult,output_energy,output_number,hits_log,count_back_whole,detector_energy_whole,hits_detectors_whole]...
                        = oneEnergyEffDistWhole(filename{i},energy_channels,outer_limit,inputfolder,detector_threshold);
                        
                    % This will be Y in our plot
                    final_Matrix2(:,i)= output_Mult;
                    
                    %final_matrix is a matrix with
                    %Energy Channel x number of different energy levels tested
                    % This will be X in our plot
                    final_Matrix(:,i) = output_energy;
                    final_Matrix3(:,i) = output_number;
                    final_Matrix4(:,i) = zeros(length(energy_channels),1);
                    final_Matrix5(i) = zeros(1,1);
                    final_Matrix6(i) = count_back_whole;
                    final_Matrix7(i) = zeros(1,1);
                    final_Matrix8(:,i) = detector_energy_whole;
                    final_Matrix9(:,i) = zeros(17,1);
                    final_Matrix10(:,i) = hits_detectors_whole;
                    final_Matrix11(:,i) = zeros(17,1);
                    
                                         
                    toc
                end
                
                %% Save Results
                %Goes to Result directory and outputs final_matrix
                cd 'E:\HERT_Drive\MATLAB Main\Result'
                save('output_MultipleParticleMatrix.mat','final_Matrix')
                disp('output_MultipleParticleMatrix.mat');
                cd ..
                
                %Goes to Efficiency_Curves Directory in prep to save
                %Eff Curve Plot
                cd Plots/Efficiency_Curves
                
                %% Obtain Data for Plotting
                r_source = 8.5;%8.5 cm for HERT-CAD
                
                % To  Set up  Plot data, we are combining the final+matrix
                %and sorting the rows
                % This is to sort the rows for plotting
                graph_data = [final_Matrix',final_Matrix2',final_Matrix3',final_Matrix4'];
                graph_data = sortrows(graph_data);
                
                graph_data2 = [final_Matrix(1,:)',final_Matrix5',final_Matrix6',final_Matrix7'];
                graph_data2 = sortrows(graph_data2);
                %After sorting, we split the matrix back to make our X and
                %Y for plotting
                graph_data_w = width(graph_data);
                %X has a column for every energy channel and rows up to the
                %number of .txt files
                x = graph_data(:,1:width(graph_data)/4);
                
                %Y has a column for every energy channel and rows up to
                %the number of .txt. files
                y_whole = graph_data(:,((width(graph_data)/4)+1):(width(graph_data)*(2/4)));
                z = graph_data(:,((width(graph_data)*2/4)+1):(3/4)*width(graph_data));
                y_inner = graph_data(:,((width(graph_data)*3/4)+1):end);
                hits_whole = y_whole;
                hits_inner = y_inner;
                
                % Sorts Out Outer Ring and Back Detector Hits
                Back_Hits_Inner = graph_data2(:,2);
                Back_Hits_Whole = graph_data2(:,3);
                Outer_Hits_Inner = graph_data2(:,4);
                
                
                
                % Calculates total number of hits across energy level and
                % energy channel for whole and inner config.
                hits_EL_whole = sum(hits_whole,2);
                hits_tot_whole = sum(hits_EL_whole);
                
                hits_EL_inner = sum(hits_inner,2);
                hits_tot_inner = sum(hits_EL_inner);
                %               part_tot_EC = z(:,:);
                
                if sim_type ==0
                    %Scales up simulated particles to total number of particles
                    part_tot_EC = 2.*z(:,:)/(1-cosd(15));
                    part_tot_EL = 2.*z(:,1)/(1-cosd(15));
                    part_tot = sum(sum(2.*z(:,:)/(1-cosd(15))));
                    
                    
                elseif sim_type ==1
                    %Full Spherical
                    part_tot_EC = z(:,:);
                    part_tot_EL = z(:,1);
                    part_tot = sum(sum(z(:,:)));
                    
                else
                    error('Error on Sim Type')
                    
                end
                %               part_tot = sum(sum(z(:,:)));
                
                %Calculates total geometric factor per energy level and
                %geometric factor of each energy channel versus incident
                %energy
                geo = (hits_tot_whole/part_tot) *(4*(pi^2)*(r_source^2));
                geo_EC = (hits_whole./part_tot_EC)*(4*(pi^2)*(r_source^2));
                geo_EL = sum(geo_EC,2);
                
                geo_inner = (hits_tot_inner/part_tot) *(4*(pi^2)*(r_source^2));
                geo_EC_inner = (hits_inner./part_tot_EC)*(4*(pi^2)*(r_source^2));
                geo_EL_inner = sum(geo_EC_inner,2);
                %geo_EL = (hits_EL./part_tot_EL)*(4*(pi^2)*(8.2^2));
                
                % Calculate Standard deviation and Error
                omega_n_whole = (part_tot_EL(:,1).*(hits_EL_whole./part_tot_EL(:,1)).*(1-(hits_EL_whole./part_tot_EL(:,1)))).^0.5;
                omega_G_whole = (4*(pi^2)*(8.2^2)).*(1-(hits_EL_whole./part_tot_EL(:,1))).*((hits_EL_whole./(part_tot_EL(:,1).^2))).^0.5;
                
                
                omega_n_inner = (part_tot_EL.*(hits_EL_inner./part_tot_EL).*(1-(hits_EL_inner./part_tot_EL(:,1)))).^0.5;
                omega_G_inner = (4*(pi^2)*(8.2^2)).*((1-(hits_EL_inner./part_tot_EL(:,1))).*(hits_EL_inner./(part_tot_EL(:,1).^2))).^0.5;
                
                
                % MeV/s Conversion Term for each detector as a function of incident energy
                graph_data3 = [final_Matrix(1,:)',final_Matrix8',final_Matrix9'];
                graph_data3 = sortrows(graph_data3);
                
                part_tot_EL_detect_whole = part_tot_EL.*ones(length(part_tot_EL),9);
                part_tot_EL_detect_inner = part_tot_EL.*ones(length(part_tot_EL),17);
                
                whole_detector_energy = graph_data3(:,2:10)./part_tot_EL_detect_whole;
                
                inner_detector_energy = graph_data3(:,11:end)./part_tot_EL_detect_inner;
                
                whole_detector_GEnergy = whole_detector_energy*(4*(pi^2)*(r_source^2));
                inner_detector_GEnergy = inner_detector_energy*(4*(pi^2)*(r_source^2));
                
                % MeV/s for each detector as a function of incident energy
                graph_data4 = [final_Matrix(1,:)',final_Matrix10',final_Matrix11'];
                graph_data4 = sortrows(graph_data4);
                
                part_tot_EL_detect_whole = part_tot_EL.*ones(length(part_tot_EL),9);
                part_tot_EL_detect_inner = part_tot_EL.*ones(length(part_tot_EL),17);
                
                whole_detector_AllCounts = graph_data4(:,2:10)./part_tot_EL_detect_whole;
                
                inner_detector_AllCounts = graph_data4(:,11:end)./part_tot_EL_detect_inner;
                
                whole_detector_GAllCounts = whole_detector_AllCounts*(4*(pi^2)*(r_source^2));
                inner_detector_GAllCounts = inner_detector_AllCounts*(4*(pi^2)*(r_source^2));
                
                % Saves variables for later graph making
                Var_String = append('OutputVariables',addin,'.mat');
                save(Var_String)
                
                addin = regexprep(addin,'_',' ');
                %% Total Geometric Factor Comparision
                line_width =2;
                textsize = 28;
                %figure
                f1 = gcf;
                f1.Position = [0 0 2000 840];
                
                hold on
                % Plot Theory Bands
                plot([x(1),x(end)],[G3_whole,G3_whole],'--g','LineWidth',line_width)
                plot([x(1),x(end)],[G3_whole_max,G3_whole_max],'--b','LineWidth',line_width)
                %plot(x(:,1),G_inner_plot,'--r','LineWidth',line_width)
                %errorbar(x(:,1),geo_EL,omega_G_whole,'b','LineWidth',line_width)
                
                %Plot Simulation Value
                plot(x(:,1),geo_EL,'-k','LineWidth',line_width)
                
                %plot(x(:,1),geo_EL,'ok','LineWidth',line_width,'MarkerFaceColor','k')
                %plot(x_FS,geo_EL_FS,'om','LineWidth',line_width,'MarkerFaceColor','m')
                
                %errorbar(x(:,1),geo_EL_inner,omega_G_inner,'m','LineWidth',line_width)
                %plot(x(:,1),geo_EL_inner,'m','LineWidth',line_width)
                
                
                %Sets yaxis to log scale. Comment out to keep plot linear
                set(gca, 'YScale', 'log')
                ylim([10^-4, 10^0])
                set(gca,'FontSize',textsize)
                %ylim([0 1])
                
                
                titlestr = append(sprintf('Total GF: %.2f MeV - %.2f MeV ',min(final_Matrix(1,:)),max(final_Matrix(1,:))),addin);
                title(titlestr,'FontSize',20)
                ylabel('Geometric Factor (cm^2 sr)','FontSize',textsize)
                xlabel('Incident Energy (MeV)','FontSize',textsize)
                
                
                %Changes legend depending on the Sim_Type
                if sim_type ==0
                    
                    legend('Theoretical Min','Theoretical Max','GEANT4 Cap','FontSize',20,'Location','southeast' )
                    
                elseif sim_type==1
                    
                    legend('Theoretical Min','Theoretical Max','GEANT4 Sphere','FontSize',20,'Location','southeast' )
                    
                end
                
                hold off
                
                %Saving the figure as a jpg then returning to main directory
                effsave = append('Total GF_',date(),'_',addin,'.jpg');
                saveas(gcf,effsave)
                
                %% Total Hits Comparision
                line_width =2;
                f2 = figure;
                f2.Position = [0 0 2000 840];
                
                hold on
                %plot(x(:,1),hits_EL_whole,'DisplayName','Whole-Simulated','LineWidth',line_width)
                %plot(x(:,1),hits_EL_inner,'DisplayName','Inner Ring-Simulated','LineWidth',line_width)
                plot(x(:,1),100*Back_Hits_Whole./output_number,'DisplayName','Whole-Back Hits Removed','LineWidth',line_width)
                %plot(x(:,1),Back_Hits_Inner,'DisplayName','Inner Ring-Back Hits Removed','LineWidth',line_width)
                %plot(x(:,1),Outer_Hits_Inner,'DisplayName','Inner Ring-Outer Ring Hits Removed','LineWidth',line_width)
                legend
                grid on
                ylim([0 100])
                yticks((0:5:100))
                titlestr = append(sprintf('Hits %.2f MeV - %.2f MeV ',min(final_Matrix(1,:)),max(final_Matrix(1,:))),addin);
                title(titlestr)
                ylabel('Percent of Hits')
                xlabel('Energy (MeV)')
                hold off
                
                %Saving the figure as a jpg then returning to main directory
                effsave = append('Total Hits',date(),addin,'.jpg');
                saveas(f2,effsave)
                
                %% Whole Ring Geometric Factor by Energy Channel
                geo_EC(geo_EC==0)=10^-31;
                line_width =2;
                textsize = 12;
                y_label = round(max(max(geo_EC)),2);
                
                if y_label < max(max(geo_EC))
                    y_label = y_label +0.01;
                end
                
                f3 = figure;
                f3.Position = [0 0 2000 840];
                EngLegend_EC = strings(1,length(energy_channels));
                color_iter = 1;
                hold on
                
                %All Channels
                channel_select = [1:length(energy_channels)];
                
                %Select which channels to highlight
                %                 channel_select = [2,10,20,30,35];
                %                 channel_select = [10];
                
                EC_plot_color = plasma(length(channel_select));
                
                hold on
                for i = 1:width(x)
                    if max(i == channel_select)
                        plot(x(:,i),geo_EC(:,i),'Color',EC_plot_color(color_iter,:),'LineWidth',line_width);
                        EngLegend_EC(i) = append(sprintf('Channel #%.0f: ',i),EngLegend(i));
                        EngLegend_EC(i) = EngLegend(i);
                        color_iter = color_iter+1;
                    else
                        plot(x(:,i),geo_EC(:,i),'Color',[0.75, 0.75, 0.75],'LineWidth',line_width);
                    end
                end
                %                 for i = 1:width(x)
                %                     if max(i == channel_select)
                %
                %                         plot(x(:,i),geo_EC(:,i),'Color',[0,0,0],'LineWidth',line_width);
                %                         %For plotting FWHM on EC plot, must calculate FWHM
                %                         %values prior to running this line
                %                         %plot([xl(i),xr(i)],[max(geo_EC(:,i))/2,max(geo_EC(:,i))/2],'-xb','LineWidth',line_width);
                %                         EngLegend_EC(i) = append(sprintf('Channel #%.0f: ',i),EngLegend(i));
                %                         EngLegend_EC(i) = EngLegend(i);
                %                         color_iter = color_iter+1;
                %                     else
                %                         plot(x(:,i),geo_EC(:,i),'Color',[0.75, 0.75, 0.75],'LineWidth',line_width);
                %                     end
                
                %end
                hold off
                
                set(gca,'FontSize',textsize)
                % plot([0,x(end)],[0.1,0.1],'--k','LineWidth',line_width)
                hold off
                titlestr_whole = append(sprintf('Geometric Factor by EC %.2f MeV - %.2f MeV ',min(final_Matrix(1,:)),max(final_Matrix(1,:))),addin);
                title(titlestr_whole,'FontSize',15)
                %title('Geometric Response Functions per Energy Channel','FontSize',15)
                ylabel('Geometric Factor (cm^2 sr)','FontSize',textsize)
                xlabel('Incident Energy (MeV)','FontSize',textsize)
                
                
                %Sets yaxis to log scale. Comment out to keep plot linear
                set(gca, 'YScale', 'log')
                
                
                ylim([10^-4, 10^0])
                grid on
                legend(EngLegend_EC,'Location','southoutside','NumColumns',8)
                
                %Saving the figure as a jpg then returning to main directory
                effsave = append('Geometric Factor Whole by EC_',date(),addin,'.jpg');
                saveas(f3,effsave)
             
                %% Efficiency
                line_width =1;
                %Create figure at a certain positon and size
                f7 = figure;
                f7.Position = [0 0 2000 840];
                hold on
                y_whole = y_whole./z;
                %Plots each energy channel with a different color
                for i = 1:width(x)
                    plot(x(:,i),y_whole(:,i),'Color',Effplotcolor(i,:),'LineWidth',line_width)
                end
                hold off
                %Adds a legend to distingish each channel
                legend(EngLegend,'Location', 'southoutside','NumColumns',8)
                
                %Adding Titles and Axis Labels
                titlestr = append(sprintf('%.2f MeV - %.2f MeV ',min(final_Matrix(1,:)),max(final_Matrix(1,:))),addin);
                title(titlestr)
                ylabel('Efficiency')
                xlabel('Energy (MeV)')
                
                %Saving the figure as a jpg then returning to main
                %directory
                effsave = append('Eff','_',addin,'.jpg');
                saveas(f7,effsave)
                cd ..
                cd ..
                
                %% Plots for Side Pen Tests
                if contains(lower(Selected_Channel_name),'side pen')
                    hits_EL_whole_1st = [0;0;6;9;8;20;19;27;50;50;78;117;132;152;180;239;280;287;344;379];
                    textsize=30;
                    sidepenFig = figure;
                    sidepenFig.Position = [0 0 1600 800];
                    hold on
                    plot(x(:,1),100*hits_EL_whole./output_number,'Color',[0 0.4470 0.7410],'LineWidth',2)
                    plot(x(:,1),100*hits_EL_whole_1st./output_number,'Color',[0.5 0 0.7410],'LineWidth',2)
                    %plot(x(:,1),100*hits_EL_inner./output_number,'Color',[0.8500 0.3250 0.0980])
                    titlestr = addin;
                    %title(titlestr,'FontSize',10)
                    ylabel('Percentage of Hits (%)','FontSize',textsize)
                    xlabel('Incident Energy (MeV)','FontSize',textsize)
                    set(gca,'FontSize',textsize)
                    grid on
                    legend("Middle Detector","1st Detector","Location","northwest",'FontSize',20)
                    
                    %Saving the figure as a jpg then returning to main
                    %directory
                    effsave = append('Side Pen',addin,'.jpg');
                    saveas(sidepenFig,effsave)
                    
                    
                    hold off
                end
                % Saves variables for later graph making
                %Var_String = append('OutputVariables',addin,'.mat');
                %save(Var_String)
                
                
                %% FWHM-Whole
                %Begins calculating FWHM values
                
                fprintf('\nWhole Configuration: Full Width at Half Max Values:\n')
                fwhm_whole = zeros(1,width(x));
                for u = 1:width(x)
                    %Using X and Y from above, calculates FWHM value for
                    %each energy channel
                    [fwhm_whole(u),xr(u),xl(u)] = findFWHM(x(:,u),geo_EC(:,u));
                    % Print full width half max values into command window
                    fprintf('%.2f - %.2f MeV: %.4f\n',energy_channels(u,1),energy_channels(u,2),fwhm_whole(u))
                end
                
                
                %Exports .txt file with the FWHM values
                cd FWHM_values
                exportname_whole = append('Whole ',inputfolder,'_',addin,'_',num2str(length(energy_channels)));
                plotsave = append(exportname_whole,'.txt');
                writematrix(fwhm_whole,plotsave)
                cd ..
                
                %Creates figure of a certain size for showing the FWHM
                %values
                b = figure;
                b.Position = [0 0 2000 840];
                hold on
                %Plots each energy channel FWHM value in the same colo as
                %the Eff. Curve
                for i = 1:size_EC(1)
                    bar((energy_channels(i,1)+energy_channels(i,2))/2,fwhm_whole(i),(energy_channels(i,2)-energy_channels(i,1)),'EdgeColor','k','FaceColor',Effplotcolor(i,:))
                    
                end
                hold off
                
                %Lables the Bar Graph
                title(titlestr_whole)
                legend(EngLegend,'Location', 'southoutside','NumColumns',8)
                ylabel('Full Width at Half Max')
                xlabel('Energy (MeV)')
                %Saves the Bar Graph as .jpg
                cd Plots/Histograms
                histsave = append(exportname_whole,'.jpg');
                saveas(gcf,histsave)
                cd ..
                cd ..
            end
    end
    choice = menu('Choose an option', 'Exit Program', 'Load one file','Load all files','Start Run');
    
end