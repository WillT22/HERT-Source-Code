%% Authored by Skyler Krantz
% Edited by Will Teague 
% Updated: Nov. 13th, 2023

% This script is the main code for converting Geant4 results to derived
% count rate and geometric factor
% Requires MatPlotLib Perceptually Uniform Colormaps

% Resets all variables and values in MATLAB
clear all;
close all;
clc;
addpath 'E:\HERT_Drive\MATLAB Main'

% Initialization of User Manipulated Variables 
detector_threshold = 0.1; % Detector Threshold (MeV)
back_threshold = 0.1; % Back Detector Threshold (MeV)
numDetect = 9;
bins = 200;
textsize = 18;

%% Calculate GEANT4 Results
% Read files from ./Result folder and store into a 1*C array
cd 'E:\HERT_Drive\MATLAB Main\Result'; % Main Result Directory

% Get Folder Names for User
topLevelFolder = pwd; % or whatever, such as 'C:\Users\John\Documents\MATLAB\work'
files = dir(topLevelFolder);
dirFlags = [files.isdir];
subFolders = files(dirFlags);
subFolderNames = {subFolders(3:end).name}; % Start at 3 to skip . and ..

FolderChoice = menu('HERT Loop: Choose an input folder', subFolderNames{:});
inputfolder = subFolderNames{FolderChoice};
addin = inputfolder;

%% Geometric Factor- Theory
% This section solves for the theoretical geometric factor of the instrument
% assuming the collimator teeth are perfect.
% Calculation derived from Sullivan 1971 paper:
% https://www.sciencedirect.com/science/article/abs/pii/0029554X71900334

% Disc 1: First Collimator Tooth
% Disc 2: Last Collimator Tooth
% Disc 3: First Detector

% Assuming the Collimator Knife Edge Stops All Particles (Min GeoFactor)
L_12 = 6.0; % cm distance between first and last collimator teeth
L_23 = 0.3; % cm distance between last collimator tooth and first detector
r1 = 0.9;   % cm radius of the first collimator tooth
r2 = 0.9;   % cm radius of the last collimator tooth
r3 = 1.0;   % cm radius of the first detector

L_13 = L_12 + L_23;

G3_whole_min = findG3whole(L_12, L_23, L_13, r1, r2, r3);

% Menu to select particle type on which max geometric factor depends
parttype_choice = menu('Particle Type', 'Electron', 'Proton');
switch parttype_choice
    case 1
        fprintf('Particle Type: Electron \n')
        % Assuming the Collimator Knife Edge Stops No Particles (Max GeoFactor)
        L_12 = 6.0; % cm distance between first and last collimator teeth
        L_23 = 0.3; % cm distance between last collimator tooth and first detector
        r1 = 1.0;   % cm radius of the first collimator tooth (larger than above)
        r2 = 1.0;   % cm radius of the last collimator tooth (larger than above)
        r3 = 1.0;   % cm radius of the first detector
        
        L_13 = L_12 + L_23;
        
        G3_whole_max = findG3whole(L_12, L_23, L_13, r1, r2, r3);
    case 2
        fprintf('Particle Type: Proton \n')
        % Protons penetrate entirety of collimator teeth
        L_12 = 6.0; % cm distance between first and last collimator teeth
        L_23 = 0.3; % cm distance between last collimator tooth and first detector
        r1 = 1.5;   % cm radius of the collimator tube interior
        r3 = 1.0;   % cm radius of the first detector
        
        L_13 = L_12 + L_23;
        
        G3_whole_max = 0.5*(pi^2)*((r1^2+r3^2+L_13^2)-(((r1^2+r3^2+L_13^2)^2-4*(r1^2)*(r3^2))^0.5));
end

%% Select Run Information
% Menu to select spherical cap or full spherical
simtype_choice = menu('Simulation Type', 'Spherical Cap (15 deg)', 'Full Sphere');
switch simtype_choice
    case 1
        sim_type = 0;
        fprintf('Sim Type: Spherical Cap (15 deg) \n')
        addin = append(addin, ' SC ');
    case 2
        sim_type = 1;
        fprintf('Sim Type: Full Sphere \n')
        addin = append(addin, ' FS ');
end

% Creates a search string for result .txt files
inputfiles = append(inputfolder, '/*.txt');
% Lists all .txt files in the Result folder
list = dir(inputfiles);
% Grabs all the names of the files in a vector (nx1 matrix)
list_fileNames = {list.name};
% Gets the number of rows and columns in file names. Columns will indicate
% the number of files that can be loaded
file_number = size(list_fileNames, 2);

% Menu to determine energy channels
channels = dir('*.txt');
energy_channel_choice = menu('Choose Energy Channels', channels.name);
energy_channels = readmatrix(channels(energy_channel_choice).name);
Selected_Channel_name = channels(energy_channel_choice).name;

% Display a menu and get a choice
choice = menu('Choose an option', 'Exit Program', 'Load one file', 'Load all files');
% Exit Program = 1 Load one file = 2 Load all files = 3 Start Run = 4
while choice ~= 1 % Choice 1 is to exit the program
    switch choice
        case 0
            disp('Error - please choose one of the options.')
        case 2 % Load One File
            % Displays all the files that can be loaded
            txt_file_choice = menu('Choose a file', list_fileNames{:});
            % If the menu button is closed, it will recycle to the initial menu
            if txt_file_choice == 0
                disp('Please select a file')
            % If the user selects a file, it will change the filename to reference that one file.
            elseif txt_file_choice > 0
                for i = 1:file_number
                    switch txt_file_choice
                        case i
                            filename = list_fileNames{i};
                    end
                end
                disp(filename); % Shows the file that was loaded
            end
        case 3 % Load all files
            filename = list_fileNames;
            fprintf('Number of files loaded: %.0f\n', length(list_fileNames))

%% Start Data Processing            
        case 4 % Start Run
            % Print detector threshold information
            detector_string = sprintf(' %.2f MeV ', detector_threshold); % Format with 2 decimal places
            addin = append(detector_string,addin);
            fprintf('Detector Threshold: %.2f MeV\n', detector_threshold)

            % Print back threshold information
            back_string = sprintf(' %.2f MeV ', back_threshold); % Format with 2 decimal places
            addin = append(back_string,addin);
            fprintf('Back Detector Threshold: %.2f MeV\n', back_threshold)
            
            % Print energy channel selection
            fprintf('%.g Energy Channels Selected \n', size(energy_channels, 1))
            addin = append(' ', erase(channels(energy_channel_choice).name, '.txt'), addin);
            
            % Creates energy channel string array for plot legend
            EngLegend = strings([1, length(energy_channels)]);
            for i = 1:size(energy_channels, 1)
                EngLegend(i) = append(num2str(round(energy_channels(i, 1), 2)), '-', num2str(round(energy_channels(i, 2), 2)), ' MeV');
            end
            
            % Creates plot colors for each energy channel
            Effplotcolor = plasma(length(energy_channels)); % Requires MatPlotLib Perceptually Uniform Colormaps
                       
            % More than one file selected
            if iscell(filename)
                disp('Start to loop oneEnergyEffDist.m');
                disp(addin);
                
                % Creates matrix to store data
                M_output_Mult = [];
                M_energy_beam = [];
                M_non_energy_beam = [];
                M_beam_number = zeros(1,file_number);
                M_back_whole = [];
                M_detector_energy_whole = zeros(9, file_number);
                M_hits_detector_whole = zeros(9, file_number);
                M_count_reject = zeros(9, file_number);
                M_run_number = zeros(1,file_number);
                
                % Nested For loops to create final matrix 1 and 2
                for i = 1:file_number
                    % For every .txt file in Results, it will run
                    % oneEnergyEffDist and add the results to finalMatrix
                    % and finalMatrix2
                    [output_Mult, energy_beam, non_energy_beam, beam_number, back_whole, hits_detectors_whole, count_reject, run_number]...
                        = oneEnergyEffDistWhole(filename{i}, energy_channels, back_threshold, inputfolder, detector_threshold);
                        
                    % This will be Y in our plot
                    M_output_Mult = [M_output_Mult, output_Mult];
                    
                    % final_matrix is a matrix with
                    % Energy Channel x number of different energy levels tested
                    % This will be X in our plot
                    M_energy_beam = [M_energy_beam, energy_beam];
                    M_non_energy_beam = [M_non_energy_beam, non_energy_beam];
                    M_beam_number(i) = beam_number;
                    M_back_whole = [M_back_whole,back_whole];
                    M_hits_detector_whole(:,i) = hits_detectors_whole;
                    M_count_reject(i) = count_reject;
                    M_run_number(i) = run_number;
                end
                
%% Obtain Data for Plotting
                % Goes to Efficiency_Curves Directory in prep to save
                % Eff Curve Plot
                cd Plots/Efficiency_Curves

                r_source = 8.5; % 8.5 cm for HERT-CAD

                % Find the indices that would sort M_output_energy
                min_incident_energy = min(M_energy_beam,[],"all");
                max_incident_energy = max(M_energy_beam,[],"all");

                % Y has a column for every energy channel and rows up to
                % the number of .txt. files
                hits_whole_EC = histcounts(M_output_Mult);
                
                 % Bin the energy for every particle simulated (hits or no)
                [M_energy_bin,M_energy_edges,M_energy_bin_indicies] = histcounts([M_energy_beam, M_non_energy_beam],bins);
                
                M_energy_midpoints = zeros(1,length(M_energy_edges)-1);
                for edge = 1:length(M_energy_edges)-1
                    M_energy_midpoints(edge) = (M_energy_edges(edge)+M_energy_edges(edge+1))/2;
                end

                % Get bin indices for all energy beam values for hit counts
                [~,~,beam_bin_indices] = histcounts(M_energy_beam, M_energy_edges);
                
                % Find the bin number for each back hit
                [~,~,back_bin_indices] = histcounts(M_back_whole, M_energy_edges);

                back_counts = zeros(1,bins);
                for bin = 1:bins
                    back_counts(bin) = nnz(back_bin_indices==bin);
                end

                if sim_type == 0
                    % Scales up simulated particles to the total number of particles
                    M_energy_bin = 2 .* M_energy_bin / (1 - cosd(15));
                    
                else
                    error('Error on Sim Type')
                end
                
                % Find the number of particles in each bin for each energy channel
                hits_log = zeros(length(energy_channels),bins);
                geo_EC = zeros(length(energy_channels),bins);
                for channel = 1:length(energy_channels)
                    for particle_index = 1:length(M_output_Mult)
                        if M_output_Mult(particle_index) == channel
                            hits_log(channel,beam_bin_indices(particle_index))= hits_log(channel,beam_bin_indices(particle_index)) + 1;
                        end
                    end
                    % Calculates the geometric factor for each channel
                    geo_EC(channel,:) = (hits_log(channel,:) ./ M_energy_bin * (4 * (pi^2) * (r_source^2)));
                end

                % Calculates total geometric factor
                geo_total = sum(geo_EC);
                
                % Saves variables for later graph making
                Var_String = append('OutputVariables', addin, '.mat');
                save(Var_String)
                
                addin = regexprep(addin, '_', ' ');
                
%% Total Geometric Factor Comparison
                line_width = 2;
                f1 = gcf;
                f1.Position = [0 0 2048 1600];
                
                hold on
                % Plot Theory Bands
                plot([10, 80], G3_whole_min * ones(1,2), '--g', 'LineWidth', line_width);
                plot([10, 80], G3_whole_max * ones(1,2), '--b', 'LineWidth', line_width);
                
                % Plot Simulation Value
                total_geo = sum(geo_EC,1);
                total_geo(total_geo < 1e-5) = 1e-5;
                
                plot(M_energy_midpoints, total_geo, '-k', 'LineWidth', line_width);

                 % Put in Penetration Limits
                xline([14,36,51],'--',{'Beryllium Window Penetration','Collimator Teeth Penetration','Veto Detector Triggering'}, ...
                    'LineWidth', 1.5,'FontSize', 19,'LabelOrientation','horizontal')
                 
                % Sets y-axis to log scale. Comment out to keep plot linear
                set(gca, 'YScale', 'log')
                xlim([10,80])
                ylim([10^-4, 10^0])
                set(gca, 'FontSize', textsize)
                
                titlestr = append(sprintf('Total GF: %.2f MeV - %.2f MeV ', min_incident_energy, max_incident_energy), addin);
                title(titlestr, 'FontSize', 20)
                ylabel('Geometric Factor (cm^2 sr)', 'FontSize', textsize)
                xlabel('Incident Energy (MeV)', 'FontSize', textsize)

                % Changes legend depending on the Sim_Type
                if sim_type == 0
                    legend_entries = {'Theoretical Min', 'Theoretical Max', 'GEANT4 Cap'};
                elseif sim_type == 1
                    legend_entries = {'Theoretical Min', 'Theoretical Max', 'GEANT4 Sphere'};
                end
                
                legend(legend_entries, 'FontSize', 20, 'Location', 'southeast');
                
                hold off
                
                % Saving the figure as a jpg then returning to the main directory
                effsave = append('Total GF_', string(datetime("today")), '_', addin, '.jpg');
                saveas(gcf, effsave)
                               
%% Total Hits Comparison
                line_width = 2;
                f2 = figure;
                f2.Position = [0 0 2048 1600];
                
                hold on
                plot(M_energy_midpoints, sum(hits_log,1) ./ M_energy_bin, 'DisplayName', 'Counted Hits', 'LineWidth', line_width)
                plot(M_energy_midpoints, back_counts ./ M_energy_bin, 'DisplayName', 'Last Detector Triggered', 'LineWidth', line_width)

                 % Put in Penetration Limits
                xline([14,35,51],'--',{'Beryllium Window Penetration','Collimator Teeth Penetration','Veto Detector Triggering'}, ...
                    'LineWidth', 1.5,'FontSize', 19,'LabelOrientation','horizontal')

                legend({'Counted Hits','Last Detector Triggered'})
                grid on
                % ylim([0 100])
                % yticks((0:5:100))
                set(gca, 'FontSize', textsize)
                titlestr = append(sprintf('Hits %.2f MeV - %.2f MeV ', min_incident_energy, max_incident_energy), addin);
                title(titlestr, 'FontSize', 20)
                ylabel('Percent of Hits')
                xlabel('Energy (MeV)')
                hold off
                
                % Saving the figure as a jpg then returning to the main directory
                effsave = append('Counted Hits', string(datetime("today")), addin, '.jpg');
                saveas(f2, effsave)

%% Geometric Factor by Energy Channel
                geo_EC(geo_EC == 0) = 10^-31;
                line_width = 2;
                y_label = round(max(max(geo_EC)), 2);
                
                if y_label < max(max(geo_EC))
                    y_label = y_label + 0.01;
                end
                
                f3 = figure;
                f3.Position = [0 0 2000 1650];
                EngLegend_EC = strings(1, length(energy_channels));
                color_iter = 1;
                hold on
                
                % All Channels
                
                % Select which channels to highlight
                % channel_select = [2,10,20,30,35];
                % channel_select = [10];
                
                hold on
                for channel = 1:size(energy_channels, 1)
                    plot(M_energy_midpoints, geo_EC(channel,:), 'Color', Effplotcolor(color_iter, :), 'LineWidth', line_width);
                    EngLegend_EC(channel) = append(sprintf('Channel #%.0f: ', channel), EngLegend(channel));
                    EngLegend_EC(channel) = EngLegend(channel);
                    color_iter = color_iter + 1;
                end
                % Put in Penetration Limits
                xline([14,35,51],'--',{'Beryllium Window Penetration','Collimator Teeth Penetration','Veto Detector Triggering'}, ...
                    'LineWidth', 1.5,'FontSize', 19,'LabelOrientation','horizontal')

                hold off
                
                set(gca, 'FontSize', textsize)
                hold off
                titlestr_whole = append(sprintf('Geometric Factor by EC %.2f MeV - %.2f MeV ', min_incident_energy, max_incident_energy), addin);
                title(titlestr_whole, 'FontSize', textsize)
                ylabel('Geometric Factor (cm^2 sr)', 'FontSize', textsize)
                xlabel('Incident Energy (MeV)', 'FontSize', textsize)
                
                % Sets y-axis to log scale. Comment out to keep the plot linear
                set(gca, 'YScale', 'log')
                
                ylim([10^-3, 10^0])
                grid on
                legend(EngLegend_EC, 'Location', 'southoutside', 'NumColumns', 6)
                
                % Saving the figure as a jpg then returning to the main directory
                effsave = append('Geometric Factor Whole by EC_', string(datetime("today")), addin, '.jpg');
                saveas(f3, effsave)
                
%% One file selected?
            else
                    disp('Start the oneEnergyEffDist.m');  
                    % Runs oneEnergyEffDistWhole for the one .txt file
                    [output_Mult,energy_beam,beam_number,back_whole,detector_energy_whole,hits_detectors_whole]...
                        = oneEnergyEffDistWhole(filename,energy_channels,back_threshold,inputfolder,detector_threshold);
                
                    hits_whole_EC = output_Mult;
                    r_source = 8.5; % 8.5 cm for HERT-CAD
                    part_tot_EC = 2 .* beam_number / (1 - cosd(15));
                    geo_EC = (hits_whole_EC ./ part_tot_EC) * (4 * (pi^2) * (r_source^2));
                    %{
                    save('output_singleParticleArray.mat','output_Mult')
                    disp('output_singleParticleArray.mat');
                    x= linspace(0,120,length(energy_channels));
                    figure
                    plot(x,output_Mult)
                    ylabel('Beam Counts')
                    xlabel('Incident Energy (MeV)')
                    %}
            end
    end 
    choice = menu('Choose an option', 'Exit Program', 'Load one file','Load all files','Start Run');
    
end
cd 'E:\HERT_Drive\MATLAB Main\Result'
close all