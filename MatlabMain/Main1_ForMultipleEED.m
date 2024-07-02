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
% Skyler had 379 descrete energy levels in his simulations (75.8 million particles)
bins = 400; % aim is to get 400 bins for comparable results
textsize = 24;
titlesize = 28;

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
L_13 = L_12 + L_23;
r1 = 0.9;   % cm radius of the first collimator tooth
r2 = 0.9;   % cm radius of the last collimator tooth
r3 = 1.0;   % cm radius of the first detector

G3_whole_min = findG3whole(L_12, L_23, L_13, r1, r2, r3);

% Menu to select particle type on which max geometric factor depends
parttype_choice = menu('Particle Type', 'Electron', 'Proton');
switch parttype_choice
    case 1
        parttype = 0;
        fprintf('Particle Type: Electron \n')
        % Assuming the Collimator Knife Edge Stops No Particles (Max GeoFactor)
        r1 = 1.0;   % cm radius of the first collimator tooth (larger than above)
        r2 = 1.0;   % cm radius of the last collimator tooth (larger than above)
        G3_whole_max = findG3whole(L_12, L_23, L_13, r1, r2, r3);

    case 2
        parttype = 1;
        fprintf('Particle Type: Proton \n')
        % Protons penetrate entirety of collimator teeth
        r1 = 1.5;   % cm radius of the collimator tube interior
        G3_whole_max = 0.5*(pi^2)*((r1^2+r3^2+L_13^2)-(((r1^2+r3^2+L_13^2)^2-4*(r1^2)*(r3^2))^0.5));
end

% r_source
r_source = 8.5; % 8.5 cm for HERT-CAD

%% Select Run Information
% Menu to select spherical cap or full spherical
simtype_choice = menu('Simulation Type', 'Spherical Cap (15 deg)', 'Full Sphere');
switch simtype_choice
    case 1
        sim_type = 0;
        fprintf('Sim Type: Spherical Cap (15 deg) \n')
        addin = append(addin, ' SC');
    case 2
        sim_type = 1;
        fprintf('Sim Type: Full Sphere \n')
        addin = append(addin, ' FS');
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
channel_path = 'channel_select';
channels = dir(fullfile(channel_path, '/*.txt'));
channel_names = {channels.name};
energy_channel_choice = menu('Choose Energy Channels', channel_names);
energy_channels = readmatrix(fullfile(channel_path, channels(energy_channel_choice).name));
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
            detector_string = sprintf(' %.2f MeV', detector_threshold); % Format with 2 decimal places
            addin = append(detector_string,addin);
            fprintf('Detector Threshold: %.2f MeV\n', detector_threshold)

            % Print back threshold information
            back_string = sprintf(' %.2f MeV', back_threshold); % Format with 2 decimal places
            addin = append(back_string,addin);
            fprintf('Back Detector Threshold: %.2f MeV\n', back_threshold)
            
            % Print energy channel selection
            fprintf(' %.g Energy Channels Selected \n', size(energy_channels, 1))
            addin = append(erase(channels(energy_channel_choice).name, '.txt'), addin);
            
            % Creates energy channel string array for plot legend
            EngLegend = strings([1, length(energy_channels)]);
            for i = 1:size(energy_channels, 1)
                EngLegend(i) = append(num2str(round(energy_channels(i, 1), 2)), '-', num2str(round(energy_channels(i, 2), 2)), ' MeV');
            end
            
            % Creates plot colors for each energy channel
            Effplotcolor = plasma(length(energy_channels)); % Requires MatPlotLib Perceptually Uniform Colormaps
            
            % Finds energy edges and midpoints based off of bins
            % logrithmic binning
                %energy_edges = logspace(log10(min_incident_energy-0.0001), log10(max_incident_energy+0.0001), bins + 1);
                %[M_energy_bin, ~ ,M_energy_bin_indicies] = histcounts([M_energy_beam, M_back_beam, M_non_energy_beam],energy_edges);
            % linear binning 
            if parttype == 0
                energy_edges = linspace(0.1,8,bins+1);
            elseif parttype == 1
                energy_edges = linspace(0,80,bins+1);
            end
            % finding midpoints
            energy_midpoints = (energy_edges(2:end) + energy_edges(1:end-1))/2;
            bin_width = energy_edges(2:end)-energy_edges(1:end-1);

            % More than one file selected
            if iscell(filename)
                disp('Start to loop oneEnergyEffDist.m');
                disp(addin);
                
                % Creates matrix to store data
                M_hit_dep = [];
                M_hit_channels = [];
                M_run_number = zeros(1,file_number);
                M_beam_number = zeros(1,file_number);
                M_energy_beam = [];
                M_non_energy_beam = [];
                M_back_beam = [];
                
                run_interest = 0;

                % Nested For loops to create final matrix 1 and 2
                for i = 1:file_number
                    % For every .txt file in Results, it will run
                    % oneEnergyEffDist and add the results to finalMatrix
                    % and finalMatrix2
                    [hit_deposited_energy, hit_energy_channels, run_number, beam_number, energy_beam, non_energy_beam, back_energy_beam]...
                        = oneEnergyEffDistWhole(filename{i}, inputfolder, energy_channels, detector_threshold, back_threshold);
                        
                    % This will be Y in our plot
                    M_hit_dep = [M_hit_dep, hit_deposited_energy'];
                    M_hit_channels = [M_hit_channels, hit_energy_channels];
                    
                    % Energy Channel x number of different energy levels tested
                    % This will be X in our plot
                    M_run_number(i) = run_number;
                    M_beam_number(i) = beam_number;
                    M_energy_beam = [M_energy_beam, energy_beam];
                    M_non_energy_beam = [M_non_energy_beam, non_energy_beam];
                    M_back_beam = [M_back_beam,back_energy_beam];

                    % To determine the file that a particle belongs
                    %{
                    if run_interest == 0 && length(M_energy_beam)>245668
                        run_interest = run_number;
                        fprintf('RUN NUMBER %.0f \n',run_number)
                    end
                    %}
                end
                
%% Obtain Data for Plotting
                % Goes to Efficiency_Curves Directory in prep to save
                % Eff Curve Plot
                cd ../Plots/Efficiency_Curves

                % Find the indices that would sort M_output_energy
                min_incident_energy = min(M_energy_beam,[],"all");
                max_incident_energy = max(M_energy_beam,[],"all");

                % Count the number of hits in each energy channel
                hits_whole_EC = histcounts(M_hit_channels,0.5:length(energy_channels)+0.5);
                
                % Bin the energy for every particle simulated (hits or no)
                [M_energy_bin,energy_edges_temp,M_energy_bin_indicies] = histcounts([M_energy_beam, M_non_energy_beam],energy_edges);
                clear energy_edges_temp;

                % Get bin indices for all energy beam values for hit counts
                [~,~,beam_bin_indices] = histcounts(M_energy_beam, energy_edges);
                
                % Find the bin number for each back hit
                [~,~,back_bin_indices] = histcounts(M_back_beam, energy_edges);

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
                hits_log = zeros(size(energy_channels,1),bins);
                geo_EC = zeros(size(energy_channels,1),bins);
                low_bins_logic = zeros(size(energy_channels,1),bins);
                
                for channel = 1:size(energy_channels,1)
                    for particle_index = 1:length(M_hit_channels)
                        if M_hit_channels(particle_index) == channel
                            hits_log(channel,beam_bin_indices(particle_index))= hits_log(channel,beam_bin_indices(particle_index)) + 1;
                        end
                    end
                    % Calculates the geometric factor for each channel
                    geo_EC(channel,:) = (hits_log(channel,:) ./ M_energy_bin * (4 * (pi^2) * (r_source^2)));

                    % Determine which bins do not contain more than two particles
                    for bin = 3:bins
                        if geo_EC(channel,bin)<10^-4 && hits_log(channel,bin-1)~=0 && hits_log(channel,bin-2)~=0
                            low_bins_logic(channel,bin) = 1;
                        end
                    end
                end
                hits_log_total = sum(hits_log,1);

                [low_bins(:,1),low_bins(:,2)] = find(low_bins_logic ~= 0);

                % Calculates total geometric factor
                geo_total = sum(geo_EC);
                
                % Saves geo_EC for later use
                %{
                fileID = fopen('geometric_factor_EC_1.txt','w');
                for channel = 1:length(energy_channels)
                    for bin = 1:bins
                        fprintf(fileID,'%.6E ',geo_EC(channel,bin));
                    end
                    fprintf(fileID,'\n');
                end
                fclose(fileID);
                %}             

                % Saves variables for later graph making
                Var_String = append('OutputVariables', addin, '.mat');
                save(Var_String)
                
                addin = regexprep(addin, '_', ' ');
                               
%% Total Hits Comparison
                line_width = 2;
                f2 = figure;
                f2.Position = [0 0 1920 1080];
                
                hold on
                plot(energy_midpoints, sum(hits_log,1) ./ M_energy_bin *100, 'DisplayName', 'Counted Hits', 'LineWidth', line_width)
                plot(energy_midpoints, back_counts ./ M_energy_bin *100, 'DisplayName', 'Last Detector Triggered', 'LineWidth', line_width)

                 % Put in Penetration Limits
                if parttype == 1
                    xline([14,35,51],'--',{'Beryllium Window Penetration','Collimator Teeth Penetration','Veto Detector Triggering'}, ...
                    'LineWidth', 1.5,'FontSize', 16,'LabelOrientation','horizontal')
                end

                legend({'Counted Hits','Last Detector Triggered'},'Location','southeast')
                grid on
                set(gca,'FontSize', textsize)
                titlestr = append(sprintf('Hits %.2f MeV - %.2f MeV ', min_incident_energy, max_incident_energy), ...
                    addin, sprintf(' %.0f Bins', bins));
                title(titlestr, 'FontSize', titlesize)
                ylabel('Percent of Hits')
                xlabel('Energy (MeV)')
                hold off
                
                % Saving the figure as a jpg then returning to the main directory
                effsave = append('Counted Hits', string(datetime("today")), addin, '.png');
                saveas(f2, effsave)

%% Total Geometric Factor Comparison
                line_width = 2;
                f1 = figure;
                f1.Position = [0 0 1800 1000];

                hold on
                % Plot Theory Bands
                plot([0,min_incident_energy,max_incident_energy], G3_whole_min * ones(1,3), '--g', 'LineWidth', line_width);
                plot([0,min_incident_energy,max_incident_energy], G3_whole_max * ones(1,3), '--b', 'LineWidth', line_width);
                
                % Plot Simulation Value
                total_geo = sum(geo_EC,1);
                total_geo(total_geo==0) = 1e-31;
                
                plot([min_incident_energy,energy_midpoints], [1e-31,total_geo], '-k', 'LineWidth', line_width);

                % Put in Penetration Limits
                %{
                if parttype == 1
                    xline([14,35,51],'--',{'Beryllium Window Penetration','Collimator Teeth Penetration','Veto Detector Triggering'}, ...
                    'LineWidth', 1.5,'FontSize', 16,'LabelOrientation','horizontal')
                end
                %}

                % Sets y-axis to log scale. Comment out to keep plot linear
                set(gca, 'YScale', 'log')
                xlim([0,max_incident_energy])
                ylim([10^-4, 10^0])
                set(gca, 'FontSize', textsize)
                
                titlestr = append(sprintf('Total GF: %.2f MeV - %.2f MeV ', min_incident_energy, max_incident_energy), ...
                    addin, sprintf(' %.0f Bins', bins));
                title(titlestr, 'FontSize', titlesize)
                %title('Proton Total Geometric Factor', 'FontSize', titlesize-2);
                ylabel('Geometric Factor (cm^2 sr)', 'FontSize', textsize)
                xlabel('Incident Energy (MeV)', 'FontSize', textsize)

                % Changes legend depending on the Sim_Type
                if sim_type == 0
                    legend_entries = {'Theoretical Min', 'Theoretical Max', 'GEANT4 Cap'};
                elseif sim_type == 1
                    legend_entries = {'Theoretical Min', 'Theoretical Max', 'GEANT4 Sphere'};
                end
                
                legend(legend_entries, 'FontSize', titlesize, 'Location', 'southeast');
                
                hold off
                
                % Saving the figure as a jpg then returning to the main directory
                effsave = append('Total GF_', string(datetime("today")), '_', addin, '.png');
                saveas(f1, effsave)

%% Geometric Factor by Energy Channel
                geo_EC(geo_EC==0)=1e-31;
                line_width = 2;
                y_label = round(max(max(geo_EC)), 2);
                
                if y_label < max(max(geo_EC))
                    y_label = y_label + 0.01;
                end
                
                f3 = figure;
                f3.Position = [0 0 1920 1080];
                hold on
                
                % All Channels
                % Select which channels to highlight
                % channel_select = [2,10,20,30,35];
                % channel_select = [1,5,10,15,18];
                % channel_select = [10];
                colors = [];
                for channel = 1:size(energy_channels, 1)
                    if exist('channel_select','var')
                        if find(channel_select == channel)>0
                            colors = [colors;plot(energy_midpoints(energy_midpoints>0.52),...
                                geo_EC(channel,energy_midpoints>0.52),...
                                'Color', Effplotcolor(channel, :),...
                                'LineWidth', line_width, 'DisplayName', EngLegend{channel})];
                        else
                            plot(energy_midpoints(energy_midpoints>0.52), geo_EC(channel,energy_midpoints>0.52), 'Color', [0.7,0.7,0.7,0.7], 'LineWidth', line_width);
                        end
                    else
                        plot(energy_midpoints(energy_midpoints>0.52), geo_EC(channel,energy_midpoints>0.52), 'Color', Effplotcolor(channel, :), 'LineWidth', line_width);
                    end
                end

                % Put in Penetration Limits
                %{
                if parttype == 1
                    xline([14,35,51],'--',{'Beryllium Window Penetration','Collimator Teeth Penetration','Veto Detector Triggering'}, ...
                    'LineWidth', 1.5,'FontSize', 16,'LabelOrientation','horizontal')
                end
                %}

                hold off
                
                set(gca, 'FontSize', textsize)
                titlestr_whole = append(sprintf('Geometric Factor by EC %.2f MeV - %.2f MeV ', min_incident_energy, max_incident_energy), ...
                    addin, sprintf(' %.0f Bins', bins));
                %title(titlestr_whole, 'FontSize', textsize)
                title('Electron Energy Channel Geometric Factor', 'FontSize', textsize)
                ylabel('Geometric Factor (cm^2 sr)', 'FontSize', textsize)
                xlabel('Incident Energy (MeV)', 'FontSize', textsize)
                
                % Sets y-axis to log scale. Comment out to keep the plot linear
                set(gca, 'YScale', 'log')
                ylim([10^-4, 10^0])
                grid on
                if exist('channel_select','var')
                    legend(colors, EngLegend(channel_select), 'Location', 'southoutside', 'NumColumns', 6);
                else
                    legend(EngLegend, 'Location', 'southoutside', 'NumColumns', 6);
                end
                
                % Saving the figure as a jpg then returning to the main directory
                effsave = append('Geometric Factor Whole by EC_', string(datetime("today")), addin, '.png');
                saveas(f3, effsave)

%% Plot counts for each energy channeldetector_threshold
%{
                f4 = figure;
                f4.Position = [0 0 1920 1080];
                hold on

                plot(1:length(energy_channels),hits_whole_EC, 'LineWidth', line_width)

                hold off
                set(gca, 'FontSize', textsize)
                ylabel('Counts', 'FontSize', textsize)
                xlabel('Energy Channel', 'FontSize', textsize)
%}

%% One file selected?
            else
                disp('Start the oneEnergyEffDist.m');  
                % Runs oneEnergyEffDistWhole for the one .txt file
                [hit_deposited_energy, hit_energy_channels, run_number, beam_number, M_energy_beam, M_non_energy_beam, back_energy_beam]...
                    = oneEnergyEffDistWhole(filename, inputfolder, energy_channels, detector_threshold, back_threshold);
                      
                hits_whole_EC = histcounts(hit_energy_channels,0.5:length(energy_channels)+0.5);
                % Get bin indices for all energy beam values for hit counts
                [~,~,beam_bin_indices] = histcounts(M_energy_beam, energy_edges);

                [M_energy_bin,energy_edges_temp,M_energy_bin_indicies] = histcounts([M_energy_beam, M_non_energy_beam],energy_edges);
                    
                if sim_type == 0
                    % Scales up simulated particles to the total number of particles
                    M_energy_bin = 2 .* M_energy_bin / (1 - cosd(15));
                    
                else
                    error('Error on Sim Type')
                end
                
                hits_log = zeros(size(energy_channels,1),bins);
                geo_EC = zeros(size(energy_channels,1),bins);
                
                for channel = 1:size(energy_channels,1)
                    for particle_index = 1:length(hit_energy_channels)
                        if hit_energy_channels(particle_index) == channel
                            hits_log(channel,beam_bin_indices(particle_index))= hits_log(channel,beam_bin_indices(particle_index)) + 1;
                        end
                    end
                    % Calculates the geometric factor for each channel
                    geo_EC(channel,:) = (hits_log(channel,:) ./ M_energy_bin * (4 * (pi^2) * (r_source^2)));
                 end
                 hits_log_total = sum(hits_log,1);

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