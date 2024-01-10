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

% Initialization of different values (I didn't know where else to put them) 
numDetect = 9;
textsize = 24;

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

% Menu to select particle type on which max geometric factor depends
parttype_choice = menu('Particle Type', 'Electron', 'Proton');
switch parttype_choice
    case 1
        sim_type = 0;
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
        sim_type = 1;
        fprintf('Particle Type: Proton \n')
        % Protons penetrate entirety of collimator teeth
        L_12 = 6.0; % cm distance between first and last collimator teeth
        L_23 = 0.3; % cm distance between last collimator tooth and first detector
        r1 = 1.5;   % cm radius of the collimator tube interior
        r3 = 1.0;   % cm radius of the first detector
        
        L_13 = L_12 + L_23;
        
        G3_whole_max = 0.5*(pi^2)*((r1^2+r3^2+L_13^2)-(((r1^2+r3^2+L_13^2)^2-4*(r1^2)*(r3^2))^0.5));
end

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
            txt_files_All = menu('Choose all files', 'All .txt files');
            switch txt_files_All
                case 1
                    % Loads all .txt files in the Result folder
                    filename = list_fileNames;
                    fprintf('Number of files loaded: %.0f\n', length(list_fileNames))
            end
        case 4 % Start Run
            % Menu to determine run type
            runtype_choice = menu('Set Back Detector Limit', '0 MeV', '0.01 MeV', '0.1 MeV');
            switch runtype_choice
                % Includes Outer Ring Hits
                case 1
                    back_limit = 0;
                    addin = append(' 0 MeV ', addin);
                case 2 % Sets 0.01 MeV Threshold for detectors
                    back_limit = 0.01;
                    addin = append(' 0.01 MeV ', addin);
                case 3 % Sets 0.1 MeV Threshold for Outer Ring
                    back_limit = 0.1;
                    addin = append(' 0.1 MeV ', addin);
            end
            % Prints Outer_limit
            fprintf('Back Detector Limit: %.2i\n', back_limit)
            
            % Menu to determine detector threshold
            detector_choice = menu('Detector Threshold', '0 MeV', '0.1 MeV');
            switch detector_choice
                case 1
                    detector_threshold = 0;
                    addin = append(addin, ' 0 MeV DT ');
                case 2
                    detector_threshold = 0.1; % MeV
                    addin = append(addin, ' 0.1 MeV DT ');
            end
            % Prints Outer_limit
            fprintf('Detector Threshold: %.2i MeV\n', detector_threshold)
            
            % Menu to determine energy channels
            channels = dir('*.txt');
            energy_channel_choice = menu('Choose Energy Channels', channels.name);
            energy_channels = readmatrix(channels(energy_channel_choice).name);
            Selected_Channel_name = channels(energy_channel_choice).name;
            fprintf('%f Energy Channels Selected \n', size(energy_channels, 1))
            addin = append(' ', erase(channels(energy_channel_choice).name, '.txt'), addin);
            
            % Preallocates EngLegend variable
            EngLegend = strings([1, length(energy_channels)]);
            % Creates Str Array for Plot legend
            for i = 1:size(energy_channels, 1)
                EngLegend(i) = append(num2str(round(energy_channels(i, 1), 2)), '-', num2str(round(energy_channels(i, 2), 2)), ' MeV');
            end
            
            % Creates n x 3 matrix for plot colors. This will give each energy channel its own color on the plot.
            % n = number of energy channels
            Effplotcolor = plasma(size(energy_channels, 1)); % Requires MatPlotLib Perceptually Uniform Colormaps
            
            % More than one file selected
            if iscell(filename)
                disp('Start to loop oneEnergyEffDist.m');
                disp(addin);
                
                % Creates matrix to store data
                M_output_energy = zeros(1, file_number);
                M_output_Mult = zeros(size(energy_channels,1), file_number);
                M_output_number = zeros(1, file_number);
                M_count_back_whole = zeros(1, file_number);
                M_detector_energy_whole = zeros(9, file_number);
                M_hits_detector_whole = zeros(9, file_number);
                
                % Nested For loops to create final matrix 1 and 2
                for i = 1: file_number
                    % For every .txt file in Results, it will run
                    % oneEnergyEffDist and add the results to finalMatrix
                    % and finalMatrix2
                    [output_Mult, output_energy, output_number, hits_log, count_back_whole, detector_energy_whole, hits_detectors_whole]...
                        = oneEnergyEffDistWhole(filename{i}, energy_channels, back_limit, inputfolder, detector_threshold);
                        
                    % This will be Y in our plot
                    M_output_Mult(:, i) = output_Mult;
                    
                    % final_matrix is a matrix with
                    % Energy Channel x number of different energy levels tested
                    % This will be X in our plot
                    M_output_energy(i) = output_energy;
                    M_output_number(i) = output_number;
                    M_count_back_whole(i) = count_back_whole;
                    M_detector_energy_whole(:, i) = detector_energy_whole;
                    M_hits_detector_whole(:, i) = hits_detectors_whole;
                end
                
                %% Save Results
                % Goes to Result directory and outputs final_matrix
                cd 'E:\HERT_Drive\MATLAB Main\Result'
                save('output_MultipleParticleMatrix.mat', 'M_output_energy')
                disp('output_MultipleParticleMatrix.mat');
                cd ..
                
                % Goes to Efficiency_Curves Directory in prep to save
                % Eff Curve Plot
                cd Plots/Efficiency_Curves
                
                %% Obtain Data for Plotting
                r_source = 8.5; % 8.5 cm for HERT-CAD

                % Find the indices that would sort M_output_energy
                [sorted_M_output_energy, sort_indices] = sort(M_output_energy);

                % Y has a column for every energy channel and rows up to
                % the number of .txt. files
                hits_whole = M_output_Mult;
                
                % Calculates total number of hits across energy level and
                % energy channel for the whole config.
                hits_EL_whole = sum(hits_whole, 1);
                hits_tot_whole = sum(hits_EL_whole);
                
                if sim_type == 0
                    % Scales up simulated particles to the total number of particles
                    part_tot_EC = 2 .* M_output_number / (1 - cosd(15));
                    part_tot_EL = part_tot_EC(:,1);
                    part_tot = sum(part_tot_EL) .* size(energy_channels, 1);
                      
                elseif sim_type == 1
                    % Full Spherical
                    part_tot_EC = M_output_number;
                    part_tot_EL = M_output_number;
                    part_tot = sum(M_output_number) .* size(energy_channels, 1);
                    
                else
                    error('Error on Sim Type')
                end
                
                % Calculates total geometric factor per energy level and
                % geometric factor of each energy channel versus incident energy
                geo = (hits_tot_whole / part_tot) * (4 * (pi^2) * (r_source^2));
                geo_EC = (hits_whole ./ (part_tot_EC .* ones(size(hits_whole)))) * (4 * (pi^2) * (r_source^2));
                geo_EL = sum(geo_EC, 1);
                
                % Calculate Standard deviation and Error
                omega_n_whole = (part_tot_EL .* (hits_EL_whole ./ part_tot_EL) .* (1 - (hits_EL_whole ./ part_tot_EL))).^0.5;
                omega_G_whole = (4 * (pi^2) * (8.2^2)) * (1 - (hits_EL_whole ./ part_tot_EL)) .* ((hits_EL_whole ./ (part_tot_EL.^2))).^0.5;
                
                % MeV/s Conversion Term for each detector as a function of incident energy (UNUSED)        
                part_tot_EL_detect_whole = part_tot_EL' .* ones(length(part_tot_EL), 9);
                whole_detector_energy = M_detector_energy_whole' ./ part_tot_EL_detect_whole;
                whole_detector_GEnergy = whole_detector_energy * (4 * (pi^2) * (r_source^2)); 
                
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
                plot([min(M_output_energy), max(M_output_energy)], G3_whole_min * ones(1,2), '--g', 'LineWidth', line_width);
                plot([min(M_output_energy), max(M_output_energy)], G3_whole_max * ones(1,2), '--b', 'LineWidth', line_width);
                
                % Plot Simulation Value
                plot(sorted_M_output_energy, geo_EL(sort_indices), '-k', 'LineWidth', line_width);

                 % Put in Penetration Limits
                xline([14,36,51],'--',{'Beryllium Window Penetration','Collimator Teeth Penetration','Veto Detector Triggering'}, ...
                    'LineWidth', 1.5,'FontSize', 19,'LabelOrientation','horizontal')
                 
                % Sets y-axis to log scale. Comment out to keep plot linear
                set(gca, 'YScale', 'log')
                ylim([10^-4, 10^0])
                set(gca, 'FontSize', textsize)
                
                titlestr = append(sprintf('Total GF: %.2f MeV - %.2f MeV ', min(M_output_energy), max(M_output_energy)), addin);
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
                plot(sorted_M_output_energy, 100 * hits_EL_whole(sort_indices) ./ M_output_number, 'DisplayName', 'Counted Hits', 'LineWidth', line_width)
                plot(sorted_M_output_energy, 100 * M_count_back_whole(sort_indices) ./ M_output_number, 'DisplayName', 'Last Detector Triggered', 'LineWidth', line_width)
                
                 % Put in Penetration Limits
                xline([14,35,51],'--',{'Beryllium Window Penetration','Collimator Teeth Penetration','Veto Detector Triggering'}, ...
                    'LineWidth', 1.5,'FontSize', 19,'LabelOrientation','horizontal')

                legend({'Counted Hits','Last Detector Triggered'})
                grid on
                % ylim([0 100])
                % yticks((0:5:100))
                set(gca, 'FontSize', textsize)
                titlestr = append(sprintf('Hits %.2f MeV - %.2f MeV ', min(M_output_energy), max(M_output_energy)), addin);
                title(titlestr, 'FontSize', 20)
                ylabel('Percent of Hits')
                xlabel('Energy (MeV)')
                hold off
                
                % Saving the figure as a jpg then returning to the main directory
                effsave = append('Counted Hits', string(datetime("today")), addin, '.jpg');
                saveas(f2, effsave)

                %% Whole Ring Geometric Factor by Energy Channel
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
                channel_select = 1:length(energy_channels);
                
                % Select which channels to highlight
                % channel_select = [2,10,20,30,35];
                % channel_select = [10];
                
                EC_plot_color = plasma(length(channel_select));
                
                hold on
                for i = 1:size(energy_channels, 1)
                    if max(i == channel_select)
                        plot(sorted_M_output_energy, geo_EC(i, sort_indices), 'Color', EC_plot_color(color_iter, :), 'LineWidth', line_width);
                        EngLegend_EC(i) = append(sprintf('Channel #%.0f: ', i), EngLegend(i));
                        EngLegend_EC(i) = EngLegend(i);
                        color_iter = color_iter + 1;
                    else
                        plot(M_output_energy(i) * ones(1, file_number), geo_EC(i, sort_indices), 'Color', [0.75, 0.75, 0.75], 'LineWidth', line_width);
                    end
                end
                % Put in Penetration Limits
                xline([14,35,51],'--',{'Beryllium Window Penetration','Collimator Teeth Penetration','Veto Detector Triggering'}, ...
                    'LineWidth', 1.5,'FontSize', 19,'LabelOrientation','horizontal')

                hold off
                
                set(gca, 'FontSize', textsize)
                hold off
                titlestr_whole = append(sprintf('Geometric Factor by EC %.2f MeV - %.2f MeV ', min(M_output_energy), max(M_output_energy)), addin);
                title(titlestr_whole, 'FontSize', textsize)
                ylabel('Geometric Factor (cm^2 sr)', 'FontSize', textsize)
                xlabel('Incident Energy (MeV)', 'FontSize', textsize)
                
                % Sets y-axis to log scale. Comment out to keep the plot linear
                set(gca, 'YScale', 'log')
                
                ylim([10^-3, 10^0])
                grid on
                legend(EngLegend_EC, 'Location', 'southoutside', 'NumColumns', 8)
                
                % Saving the figure as a jpg then returning to the main directory
                effsave = append('Geometric Factor Whole by EC_', string(datetime("today")), addin, '.jpg');
                saveas(f3, effsave)
                
%{
                %% Efficiency
                line_width = 1;
                % Create figure at a certain position and size
                f7 = figure;
                f7.Position = [0 0 2000 840];
                hold on
                hits_rate = hits_whole ./ M_output_number;
                % Plots each energy channel with a different color
                for i = 1:size(energy_channels, 1)
                    plot(sorted_M_output_energy, hits_rate(i, sort_indices), 'Color', Effplotcolor(i, :), 'LineWidth', line_width)
                end
                hold off
                
                % Adds a legend to distinguish each channel
                legend(EngLegend, 'Location', 'southoutside', 'NumColumns', 8)
                
                % Adding Titles and Axis Labels
                titlestr = append(sprintf('%.2f MeV - %.2f MeV ', min(M_output_energy), max(M_output_energy)), addin);
                title(titlestr)
                ylabel('Efficiency')
                xlabel('Energy (MeV)')
                
                % Saving the figure as a jpg then returning to the main directory
                effsave = append('Eff', '_', addin, '.jpg');
                saveas(f7, effsave)
                cd ..
                cd ..
%}
%{
                %% Plots for Side Pen Tests
                % Check if the selected energy channel name contains 'side pen'
                if contains(lower(Selected_Channel_name),'side pen')
                    % Define simulated hits for the first detector in side pen tests
                    hits_EL_whole_1st = [0;0;6;9;8;20;19;27;50;50;78;117;132;152;180;239;280;287;344;379];
                    
                    % Set figure properties
                    textsize=30;
                    sidepenFig = figure;
                    sidepenFig.Position = [0 0 1600 800];
                    
                    % Plot simulated hits for the middle and first detector
                    hold on
                    plot(M_output_energy,100*hits_EL_whole./M_output_number,'Color',[0 0.4470 0.7410],'LineWidth',2)
                    plot(M_output_energy,100*hits_EL_whole_1st./M_output_number,'Color',[0.5 0 0.7410],'LineWidth',2)
                    
                    % Set plot properties
                    titlestr = addin;
                    ylabel('Percentage of Hits (%)','FontSize',textsize)
                    xlabel('Incident Energy (MeV)','FontSize',textsize)
                    set(gca,'FontSize',textsize)
                    grid on
                    legend("Middle Detector","1st Detector","Location","northwest",'FontSize',20)
                    
                    % Save the figure as a jpg
                    effsave = append('Side Pen',addin,'.jpg');
                    saveas(sidepenFig,effsave)
                    
                    % Release figure hold
                    hold off
                end
%}                
%{
                %% FWHM-Whole
                % Calculate Full Width at Half Max (FWHM) values for each energy channel
                fprintf('\nWhole Configuration: Full Width at Half Max Values:\n')
                fwhm_whole = zeros(size(energy_channels));
                for u = 1:size(energy_channels,1)
                    % Using X and Y from above, calculates FWHM value for each energy channel
                    [fwhm_whole(u),xr(u),xl(u)] = findFWHM(M_output_energy,geo_EC(:,u));
                    % Print full width half max values into the command window
                    fprintf('%.2f - %.2f MeV: %.4f\n',energy_channels(u,1),energy_channels(u,2),fwhm_whole(u))
                end
                
                % Exports .txt file with the FWHM values
                cd FWHM_values
                exportname_whole = append('Whole ',inputfolder,'_',addin,'_',num2str(length(energy_channels)));
                plotsave = append(exportname_whole,'.txt');
                writematrix(fwhm_whole,plotsave)
                cd ..
                
                % Create figure for showing the FWHM values
                b = figure;
                b.Position = [0 0 2000 840];
                hold on
                
                % Plot each energy channel's FWHM value in the same color as the Eff. Curve
                for i = 1:size(energy_channels,1)
                    bar((energy_channels(i,1)+energy_channels(i,2))/2,fwhm_whole(i),(energy_channels(i,2)-energy_channels(i,1)),'EdgeColor','k','FaceColor',Effplotcolor(i,:))
                end
                hold off
                
                % Label the Bar Graph
                title(titlestr_whole)
                legend(EngLegend,'Location', 'southoutside','NumColumns',8)
                ylabel('Full Width at Half Max')
                xlabel('Energy (MeV)')
                
                % Save the Bar Graph as .jpg
                cd Plots/Histograms
                histsave = append(exportname_whole,'.jpg');
                saveas(gcf,histsave)
                cd ..
                cd ..
%} 

            %Diagnostics for one file selected
            else
                    disp('Start the oneEnergyEffDist.m');  
                    % Runs oneEnergyEffDistWhole for the one .txt file
                    [output_Mult,output_energy,output_number,hits_log,count_back_whole,detector_energy_whole,hits_detectors_whole]...
                        = oneEnergyEffDistWhole(filename,energy_channels,back_limit,inputfolder,detector_threshold);
                
                    hits_whole = output_Mult;
                    r_source = 8.5; % 8.5 cm for HERT-CAD
                    part_tot_EC = 2 .* output_number / (1 - cosd(15));
                    geo_EC = (hits_whole ./ part_tot_EC) * (4 * (pi^2) * (r_source^2));
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
close all