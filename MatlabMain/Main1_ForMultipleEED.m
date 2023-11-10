%%HERT: mainForMultipleEED.mfiletype
% Last modified: 11/9/2023


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

% Disc 1: First Collimator Tooth
% Disc 2: Last Collimator Tooth
% Disc 3: First Detector

% Assuming the Collimator Knife Edge Stops All Particles (Min GeoFactor)
L_12 = 6.0;%cm distance between first and last collimator teeth
L_23 = 0.3; %cm distance between last collimator tooth and first detector
r1 = 0.9; %cm radius of first collimator tooth
r2 = 0.9; %cm radius of last collimator tooth
r3 = 1.0; %cm radius of first detector

L_13 = L_12+L_23;

G3_whole_min = findG3whole(L_12,L_23,L_13,r1,r2,r3);

% Assuming the Collimator Knife Edge Stops No Particles (Max GeoFactor)
L_12 = 6.0;%cm distance between first and last collimator teeth
L_23 = 0.3; %cm distance between last collimator tooth and first detector
r1 = 1.0; %cm radius of first collimator tooth (larger than above)
r2 = 1.0; %cm radius of last collimator tooth (larger than above)
r3 = 1.0; %cm radius of first detector

L_13 = L_12+L_23;

G3_whole_max = findG3whole(L_12,L_23,L_13,r1,r2,r3);

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
file_number = size(list_fileNames,2);

% Display a menu and get a choice
choice = menu('Choose an option', 'Exit Program', 'Load one file','Load all files');
%Exit Program =1 Load one file = 2 Load all files =3 Start Run=4
while choice ~= 1 % Choice 1 is to exit the program
    switch choice
        case 0
            disp('Error - please choose one of the options.')
        case 2 %Load One File
            %Displays all the files that can be loaded
            txt_file_choice = menu('Choose a file',list_fileNames{:});
            %If menu button is closed, it will recycle to inital menu
            if txt_file_choice == 0
                disp('Please select a file')
            %If user selects a file, it will change filename to reference that one file.
            elseif txt_file_choice > 0
                for i = 1:file_number
                    switch txt_file_choice
                        case i
                            filename = list_fileNames{i};
                    end
                end
                disp(filename); %Shows the file that was loaded
            end
        case 3 %Load all files
            txt_files_All= menu('Choose all files','All .txt files');
            switch txt_files_All
                case 1
                    %Loads all .txt files in Result folder
                    filename = list_fileNames';
                    fprintf('Number of files loaded: %.0f\n',length(list_fileNames))       
            end
        case 4 % Start Run
            %Menu to determine runtype
            runtype_choice = menu('Set Back Detector Limit','0 MeV','0.01 MeV','0.1 MeV');
            switch runtype_choice
                %Includes Outer Ring Hits
                case 1
                    back_limit = 0;
                    addin = append(' 0 MeV ',addin);
                case 2 %Sets 0.01 MeV Threshhold for detectors
                    back_limit = 0.01;
                    addin = append(' 0.01 MeV ',addin);
                case 3 %Sets 0.1 MeV Threshhold for Outer Ring
                    back_limit = 0.1;
                    addin = append(' 0.1 MeV ',addin);        
            end
            %Prints Outer_limit
            fprintf('Back Detector Limit: %.2i\n',back_limit)
            
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
            fprintf('Detector Threshold: %.2i MeV\n',detector_threshold)
            
            %Menu to determine energy channels
            %Finds all .txt files in Results directory. These txt files contain the energy channels
            channels = dir('*.txt');
            
            %Produces menu for user to select while energy channels
            energy_channel_choice = menu('Choose Energy Channels',channels.name);
            
            %Reads in selected energy channel file and prepares addin for result plots
            energy_channels = readmatrix(channels(energy_channel_choice).name);
            Selected_Channel_name = channels(energy_channel_choice).name;
            fprintf('%f Energy Channels Selected \n',size(energy_channels,1))
            addin = append(' ',erase(channels(energy_channel_choice).name,'.txt'),addin);
            
            %Preallocates EngLegend variable
            EngLegend = strings([1,length(energy_channels)]);
            %Creates Str Array for Plot legend
            for i = 1:size(energy_channels,1)
                EngLegend(i) = append(num2str(round(energy_channels(i,1),2)),'-',num2str(round(energy_channels(i,2),2)),' MeV');
            end
            
            %Creates n x 3 matrix for plot colors. This will give each energy channel its own color on the plot.
            %n = number of energy channels
            Effplotcolor = plasma(size(energy_channels,1)); %requires MatPlotLib Perceptually Uniform Colormaps

%%THIS SECTION DOES NOT WORK, WHAT IS OUTPUT0           
            %One file selected
            if size(filename,1)==1
                disp('Start the oneEnergyEffDist.m');  
                %Runs oneEnergyEffDistWhole for the one .txt file
                [output_Mult,output_energy,output_number,hits_log,count_back_whole,detector_energy_whole,hits_detectors_whole]...
                    = oneEnergyEffDistWhole(filename,energy_channels,back_limit,inputfolder,detector_threshold);
                    
                hits_whole = sum(output0.*hits_whole);
                save('output_singleParticleArray.mat','output0')
                disp('output_singleParticleArray.mat');
                x= linspace(1.0,7.0,length(energy_channels));
                hits_whole = output0';
                figure
                plot(x,hits_whole(:,:))

%%THIS SECTION STARTS WORKING AGAIN                
            % More than one file selected
            elseif size(filename,1) > 1
                disp('Start to loop oneEnergyEffDist.m');
                disp(addin);
                
                %Creates matrix to store data
                M_output_energy = zeros(file_number,1);
                M_output_Mult = zeros(length(energy_channels),file_number);
                M_output_number = zeros(length(energy_channels),file_number);
                M_count_back_whole = zeros(1,file_number);
                M_detector_energy_whole = zeros(9,file_number);
                M_hits_detector_whole = zeros(9,file_number);
                
                %Nested For loops to create final matrix 1 and 2
                for i = 1: file_number
                    %For every .txt file in Results, it will run
                    %oneEnergyEffDist and add the results to finalMatrix
                    %and finalMatrix2
                    [output_Mult,output_energy,output_number,hits_log,count_back_whole,detector_energy_whole,hits_detectors_whole]...
                        = oneEnergyEffDistWhole(filename{i},energy_channels,back_limit,inputfolder,detector_threshold);
                        
                    % This will be Y in our plot
                    M_output_Mult(:,i)= output_Mult;
                    
                    %final_matrix is a matrix with
                    %Energy Channel x number of different energy levels tested
                    % This will be X in our plot
                    M_output_energy(i) = output_energy;
                    M_output_number(:,i) = output_number;
                    M_count_back_whole(i) = count_back_whole;
                    M_detector_energy_whole(:,i) = detector_energy_whole;
                    M_hits_detector_whole(:,i) = hits_detectors_whole;
                end
                
                %% Save Results
                %Goes to Result directory and outputs final_matrix
                cd 'E:\HERT_Drive\MATLAB Main\Result'
                save('output_MultipleParticleMatrix.mat','M_output_energy')
                disp('output_MultipleParticleMatrix.mat');
                cd ..
                
                %Goes to Efficiency_Curves Directory in prep to save
                %Eff Curve Plot
                cd Plots/Efficiency_Curves
                
                %% Obtain Data for Plotting
                r_source = 8.5; %8.5 cm for HERT-CAD

                %Y has a column for every energy channel and rows up to
                %the number of .txt. files
                hits_whole = M_output_Mult';
               
                % Sorts Out Outer Ring and Back Detector Hits
                Back_Hits_Whole = M_count_back_whole';
                
                % Calculates total number of hits across energy level and
                % energy channel for whole config.
                hits_EL_whole = sum(hits_whole,2);
                hits_tot_whole = sum(hits_EL_whole);
                
                if sim_type ==0
                    %Scales up simulated particles to total number of particles
                    part_tot_EC = 2.*M_output_number'/(1-cosd(15));
                    part_tot_EL = 2.*M_output_number(1,:)'/(1-cosd(15));
                    part_tot = sum(sum(2.*M_output_number'/(1-cosd(15))));
                      
                elseif sim_type ==1
                    %Full Spherical
                    part_tot_EC = M_output_number';
                    part_tot_EL = M_output_number(1,:)';
                    part_tot = sum(sum(M_output_number,2));
                    
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
                
                % Calculate Standard deviation and Error
                omega_n_whole = (part_tot_EL.*(hits_EL_whole./part_tot_EL).*(1-(hits_EL_whole./part_tot_EL))).^0.5;
                omega_G_whole = (4*(pi^2)*(8.2^2)).*(1-(hits_EL_whole./part_tot_EL)).*((hits_EL_whole./(part_tot_EL.^2))).^0.5;
                
                % MeV/s Conversion Term for each detector as a function of incident energy        
                part_tot_EL_detect_whole = part_tot_EL.*ones(length(part_tot_EL),9);
                
                whole_detector_energy = M_detector_energy_whole'./part_tot_EL_detect_whole;
                
                whole_detector_GEnergy = whole_detector_energy*(4*(pi^2)*(r_source^2)); 
                
                % MeV/s for each detector as a function of incident energy
                part_tot_EL_detect_whole = part_tot_EL.*ones(length(part_tot_EL),9);
                
                whole_detector_AllCounts = M_hits_detector_whole'./part_tot_EL_detect_whole;
                
                whole_detector_GAllCounts = whole_detector_AllCounts*(4*(pi^2)*(r_source^2));
                
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
                plot(M_output_energy,G3_whole_min*ones(length(M_output_energy)),'--g','LineWidth',line_width)
                plot(M_output_energy,G3_whole_max*ones(length(M_output_energy)),'--b','LineWidth',line_width)

                %Plot Simulation Value
                plot(M_output_energy,geo_EL,'-k','LineWidth',line_width)
                 
                %Sets yaxis to log scale. Comment out to keep plot linear
                set(gca, 'YScale', 'log')
                ylim([10^-4, 10^0])
                set(gca,'FontSize',textsize)
                %ylim([0 1])
                
                
                titlestr = append(sprintf('Total GF: %.2f MeV - %.2f MeV ',min(M_output_energy),max(M_output_energy)),addin);
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
                effsave = append('Total GF_',string(datetime("today")),'_',addin,'.jpg');
                saveas(gcf,effsave)
                
                %% Total Hits Comparision
                line_width =2;
                f2 = figure;
                f2.Position = [0 0 2000 840];
                
                hold on
                plot(M_output_energy,100*Back_Hits_Whole./output_number,'DisplayName','Whole-Back Hits Removed','LineWidth',line_width)
                legend
                grid on
                ylim([0 100])
                yticks((0:5:100))
                titlestr = append(sprintf('Hits %.2f MeV - %.2f MeV ',min(M_output_energy),max(M_output_energy)),addin);
                title(titlestr)
                ylabel('Percent of Hits')
                xlabel('Energy (MeV)')
                hold off
                
                %Saving the figure as a jpg then returning to main directory
                effsave = append('Total Hits',string(datetime("today")),addin,'.jpg');
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
                channel_select = 1:length(energy_channels);
                
                %Select which channels to highlight
                %                 channel_select = [2,10,20,30,35];
                %                 channel_select = [10];
                
                EC_plot_color = plasma(length(channel_select));
                
                hold on
                for i = 1:size(energy_channels,1)
                    if max(i == channel_select)
                        plot(M_output_energy,geo_EC(:,i),'Color',EC_plot_color(color_iter,:),'LineWidth',line_width);
                        EngLegend_EC(i) = append(sprintf('Channel #%.0f: ',i),EngLegend(i));
                        EngLegend_EC(i) = EngLegend(i);
                        color_iter = color_iter+1;
                    else
                        plot(M_output_energy(i)*ones(1,file_number),geo_EC(:,i),'Color',[0.75, 0.75, 0.75],'LineWidth',line_width);
                    end
                end
                hold off
                
                set(gca,'FontSize',textsize)
                hold off
                titlestr_whole = append(sprintf('Geometric Factor by EC %.2f MeV - %.2f MeV ',min(M_output_energy),max(M_output_energy)),addin);
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
                effsave = append('Geometric Factor Whole by EC_',string(datetime("today")),addin,'.jpg');
                saveas(f3,effsave)
             
                %% Efficiency
                line_width =1;
                %Create figure at a certain positon and size
                f7 = figure;
                f7.Position = [0 0 2000 840];
                hold on
                hits_whole = hits_whole./(M_output_number');
                %Plots each energy channel with a different color
                for i = 1:size(energy_channels,1)
                    plot(M_output_energy,hits_whole(:,i),'Color',Effplotcolor(i,:),'LineWidth',line_width)
                end
                hold off
                %Adds a legend to distingish each channel
                legend(EngLegend,'Location', 'southoutside','NumColumns',8)
                
                %Adding Titles and Axis Labels
                titlestr = append(sprintf('%.2f MeV - %.2f MeV ',min(M_output_energy),max(M_output_energy)),addin);
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
                    plot(M_output_energy,100*hits_EL_whole./output_number,'Color',[0 0.4470 0.7410],'LineWidth',2)
                    plot(M_output_energy,100*hits_EL_whole_1st./output_number,'Color',[0.5 0 0.7410],'LineWidth',2)
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
                fwhm_whole = zeros(size(energy_channels));
                for u = 1:size(energy_channels,1)
                    %Using X and Y from above, calculates FWHM value for
                    %each energy channel
                    [fwhm_whole(u),xr(u),xl(u)] = findFWHM(M_output_energy,geo_EC(:,u));
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
                for i = 1:size(energy_channels,1)
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
close all