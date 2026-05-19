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
addpath 'D:\HERT_Drive\MATLAB Main'

% Initialization of User Manipulated Variables 
detector_threshold = 0.1; % Detector Threshold (MeV)
numDetect = 9;
textsize = 20;
titlesize = 28;

%% Calculate GEANT4 Results
% Read files from ./Result folder and store into a 1*C array
cd 'D:\HERT_Drive\MATLAB Main\Result'; % Main Result Directory

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
        back_threshold = 0.1; % Back Detector Threshold (MeV)
    case 2
        parttype = 1;
        fprintf('Particle Type: Proton \n')
        % Protons penetrate entirety of collimator teeth
        r1 = 1.5;   % cm radius of the collimator tube interior
        G3_whole_max = 0.5*(pi^2)*((r1^2+r3^2+L_13^2)-(((r1^2+r3^2+L_13^2)^2-4*(r1^2)*(r3^2))^0.5));
        
        r1 = 1; % cm radius of the first detector
        r9 = 2; % cm radius of the last detector
        L_19 = 0.15 * 7 + 0.108 * 4 + 0.1 * 4; % distance between back of first detector and front of last detector
        G3_pen_theory = 0.5*(pi^2)*((r1^2+r9^2+L_19^2)-(((r1^2+r9^2+L_19^2)^2-4*(r1^2)*(r9^2))^0.5)) * 2;
        back_threshold = 0.1; % Back Detector Threshold (MeV)
end

% r_source
r_source = 8.5; % 8.5 cm for HERT-CAD
% r_source = 6.45736; % 6.45736 cm for HERT_DARTBe
source_angle = 15; % 15 for spherical cap, 20.825/2 for HERT testing

%% Select Run Information
% Menu to select spherical cap or full spherical
source_label = sprintf('Spherical Cap (%.0f deg)', source_angle);
simtype_choice = menu('Simulation Type', source_label, 'Full Sphere');
switch simtype_choice
    case 1
        sim_type = 0;
        fprintf('Sim Type: %d.\n', simtype_choice)
        addin = append(addin, ' SC');
    case 2
        sim_type = 1;
        fprintf('Sim Type: %d.\n', simtype_choice)
        addin = append(addin, ' FS');
end

% Creates a search string for result .txt files
file_prefix = 'PostProcessHERT';
file_suffix = '.txt';
inputfiles = append(inputfolder, '/', file_prefix, '*', file_suffix);
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
fprintf('%.0f Energy Channels Selected \n', size(energy_channels, 1))
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
            addin = append(erase(channels(energy_channel_choice).name, '.txt'), addin);
                
            % Finds energy edges and midpoints based off of bins
            if parttype == 0
                bins = 400;
                energy_min = 0.01;
                energy_max = 10;
                energy_edges = logspace(log10(energy_min),log10(energy_max),bins+1);
                energy_edges(end) = energy_max + 0.01;
            elseif parttype == 1
                bins = 400;
                energy_min = 10;
                energy_max = 1000;
                energy_edges = logspace(log10(energy_min),log10(energy_max),bins+1);
            end
            % finding midpoints
            energy_midpoints = (energy_edges(2:end) + energy_edges(1:end-1))/2;
            bin_width = diff(energy_edges);

            % More than one file selected
            if iscell(filename)
                disp('Start to loop oneEnergyEffDist.m');
                disp(addin);
                
                % Creates matrix to store data
                M_hit_dep = [];
                M_hit_channels = [];
                M_hit_detectors = zeros(1,9);

                M_back_dep = [];
                M_back_channels = [];
                M_back_detectors = zeros(1,9);

                M_run_number = zeros(1,file_number);
                M_beam_number = zeros(1,file_number);
                M_energy_beam = [];
                M_non_energy_beam = [];
                M_run_number_back = zeros(1,file_number);
                M_beam_number_back = zeros(1,file_number);
                M_back_beam = [];
                M_energy_bin = zeros(1,bins);
                min_incident_energy = 1000;
                max_incident_energy = 0;
                run_interest = [];
                run_interest_back = [];

                % Nested For loops to create final matrix 1 and 2
                for i = 1:file_number
                    % For every .txt file in Results, it will run
                    % oneEnergyEffDist and add the results to finalMatrix
                    % and finalMatrix2
                    [hit_deposited_energy, hit_energy_channels, hit_detectors, energy_beam, run_number, beam_number, ...
                        back_deposited_energy, back_energy_channels, back_detectors, back_energy_beam, non_energy_beam]...
                        = oneEnergyEffDistWhole(filename{i}, inputfolder, energy_channels, detector_threshold, back_threshold);
                        
                    % This will be Y in our plot
                    M_hit_dep = [M_hit_dep, hit_deposited_energy'];
                    M_hit_channels = [M_hit_channels, hit_energy_channels];
                    M_hit_detectors = M_hit_detectors + hit_detectors;

                    M_back_dep = [M_back_dep, back_deposited_energy'];
                    M_back_channels = [M_back_channels, back_energy_channels];
                    M_back_detectors = M_back_detectors + back_detectors;
                    
                    % Energy Channel x number of different energy levels tested
                    % This will be X in our plot
                    M_energy_beam = [M_energy_beam, energy_beam];
                    M_run_number(i) = run_number;
                    M_beam_number(i) = beam_number;
                    run_interest = [run_interest, ones(1,length(energy_beam))*run_number];
                    % This one takes a lot of memory
                    %M_non_energy_beam = [M_non_energy_beam, non_energy_beam];
                    M_back_beam = [M_back_beam,back_energy_beam];
                    M_run_number_back(i) = run_number;
                    M_beam_number_back(i) = beam_number;
                    run_interest_back = [run_interest_back, ones(1,length(back_energy_beam))*run_number];

                    % Bin the energy for every particle simulated (hits or no)
                    [M_energy_bin_temp,~,M_energy_bin_indicies] = histcounts([energy_beam, non_energy_beam],energy_edges);
                    M_energy_bin = M_energy_bin + M_energy_bin_temp;

                    % Find the indices that would sort M_output_energy
                    if min_incident_energy > min(energy_beam,[],"all");
                        min_incident_energy = min(energy_beam,[],"all");
                    end
                    if max_incident_energy < max(energy_beam,[],"all")
                        max_incident_energy = max(energy_beam,[],"all");
                    end
                    
                    clear energy_beam
                    clear non_energy_beam
                    clear back_energy_beam
                end

                particles = [M_energy_beam',run_interest'];
                back_particles = [M_back_beam',run_interest_back'];

%% Obtain Data for Plotting
                % Goes to Efficiency_Curves Directory in prep to save
                % Eff Curve Plot
                cd ../Plots/Efficiency_Curves

                % Count the number of hits in each energy channel
                hits_whole_EC = histcounts(M_hit_channels,0.5:length(energy_channels)+0.5);
                back_whole_EC = histcounts(M_back_channels,0.5:length(energy_channels)+0.5);
 
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
                    M_energy_bin = 2 .* M_energy_bin / (1 - cosd(source_angle));
                    M_back_beam  = 2 .* M_back_beam  / (1 - cosd(source_angle));
                elseif sim_type == 1
                else
                    error('Error on Sim Type')
                end

                % Find the number of particles in each bin for each energy channel
                hits_log = zeros(size(energy_channels,1),bins);
                back_log = zeros(size(energy_channels,1),bins);
                geo_EC = zeros(size(energy_channels,1),bins);
                geo_back_EC = zeros(size(energy_channels,1),bins);
                low_bins_logic = zeros(size(energy_channels,1),bins);
                
                for channel = 1:size(energy_channels,1)
                    for particle_index = 1:length(M_hit_channels)
                        if M_hit_channels(particle_index) == channel
                            hits_log(channel,beam_bin_indices(particle_index))= hits_log(channel,beam_bin_indices(particle_index)) + 1;
                        end
                    end
                    for particle_index_back = 1:length(M_back_channels)
                        if M_back_channels(particle_index_back) == channel
                            back_log(channel,back_bin_indices(particle_index_back))= back_log(channel,back_bin_indices(particle_index_back)) + 1;
                        end
                    end
                    % Calculates the geometric factor for each channel
                    geo_EC(channel,:) = hits_log(channel,:) ./ M_energy_bin * (4 * (pi^2) * (r_source^2));
                    geo_back_EC(channel,:) = back_log(channel,:) ./ M_energy_bin * (4 * (pi^2) * (r_source^2));
                end
                if parttype == 1
                    geo_back_EC_combined = [sum(geo_back_EC(1:8, :), 1); geo_back_EC(9:end, :)];
                    energy_channels_combined = [energy_channels(1, 1), energy_channels(8, 2); energy_channels(9:end, :)];
                end
                
                hits_log_total = sum(hits_log,1);
                back_log_total = sum(back_log,1);

                [low_bins(:,1),low_bins(:,2)] = find(low_bins_logic ~= 0);

                % Saves variables for later graph making
                Var_String = append('OutputVariables', addin, '.mat');
                save(Var_String)
                
                addin = regexprep(addin, '_', ' ');
                               
%% Total Hits Comparison
%{
                line_width = 2;
                f2 = figure;
                f2.Position = [0 0 1920 1080];
                
                hold on
                plot(energy_midpoints, hits_log_total ./ M_energy_bin *100, 'DisplayName', 'Counted Hits', 'LineWidth', line_width)
                plot(energy_midpoints, back_counts ./ M_energy_bin * 100, 'DisplayName', 'Last Detector Triggered', 'LineWidth', line_width)

                % Put in Penetration Limits
                if parttype == 1
                    xline([14.15,52.29,69.09],'--',{'Be Penetration','D9 Triggering', 'W Penetration'}, ...
                    'LineWidth', 1.5,'FontSize', 16,'LabelOrientation','aligned')
                end

                legend({'Counted Hits','Last Detector Triggered'},'Location','northeast')
                grid on
                set(gca,'FontSize', textsize)
                titlestr = append(sprintf('Hits %.2f MeV - %.2f MeV ', min_incident_energy, max_incident_energy), ...
                    addin, sprintf(' %.0f Bins', bins));
                %title(titlestr, 'FontSize', titlesize)
                if parttype == 0
                    set(gca, 'XTick', [0.5,1,1.5,2,3,4,5,6,7,8]);
                    xlim([0.4,max_incident_energy])
                    ylim([10^-5, 10^-1])
                elseif parttype == 1
                    set(gca, 'XTick', [0,10,20,30,40,50,60,80,100,120,140,160,200,300,400,500,600,800,1000]);
                    xlim([10,max_incident_energy])
                    ylim([10^-3, 10^1])
                end
                set(gca, 'XScale', 'log')
                set(gca, 'YScale', 'log')
                ylabel('Percent of Hits')
                xlabel('Energy (MeV)')
                hold off
                
                % Saving the figure as a jpg then returning to the main directory
                effsave = append('Counted Hits', string(datetime("today")), addin, '.png');
                saveas(f2, effsave)
%}

%% Total Geometric Factor Comparison
back_bool = exist('geo_back_EC_combined','var');

line_width = 2;
f1 = figure;
f1.Position = [0 0 1800 1000];

hold on
% Plot Theory Bands
%plot([0,min_incident_energy,max_incident_energy], G3_whole_min * ones(1,3), '--g', 'LineWidth', line_width);
%plot([0,min_incident_energy,max_incident_energy], G3_whole_max * ones(1,3), '--b', 'LineWidth', line_width);

% Plot Simulation Value
total_geo = sum(geo_EC,1);
total_geo(total_geo==0) = 1e-31;

plot([min_incident_energy,energy_midpoints], [1e-31,total_geo], 'LineWidth', line_width);

if back_bool
    back_total_geo = sum(geo_back_EC,1);
    back_total_geo(back_total_geo==0) = 1e-31;

    plot([min_incident_energy,energy_midpoints], [1e-31,back_total_geo], 'LineWidth', line_width);
end

%plot([min_incident_energy,energy_midpoints_noveto], [1e-31,total_geo_noveto], '-r', 'LineWidth', line_width);

% Put in Penetration Limits
if parttype == 0
    xline([0.62, 5.71, 15.91, 30.89], '--', {'Be Penetration','D9 Penetration', 'W Penetration', 'W & Stack Penetration'}, ...
        'LineWidth', 1.5, 'FontSize', 16, 'LabelOrientation', 'aligned');
elseif parttype == 1
    xline([14.15,52.29,69.09,99.64],'--',{'Be Penetration','D9 Triggering', 'W Penetration', 'W & Stack Penetration'}, ...
    'LineWidth', 1.5,'FontSize', 16,'LabelOrientation','aligned')

    % Plot theory band
    yline([G3_pen_theory],'--',{'Proton Penetration Theory'}, ...
    'LineWidth', 1.5,'FontSize', 16,'LabelOrientation','aligned','LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment','bottom')
end

legend({'Range','Penetrating'},'Location','northeast')
grid on
% Sets y-axis to log scale. Comment out to keep plot linear
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
if parttype == 0
    set(gca, 'XTick', [0.5,1,1.5,2,3,4,5,6,7,8]);
    xlim([0.4,max_incident_energy])
    ylim([10^-3, 10^0])
elseif parttype == 1
    set(gca, 'XTick', [0,10,20,30,40,50,60,80,100,120,150,200,300,400,500,600,800,1000]);
    xlim([10,max_incident_energy])
    ylim([10^-2, 10^2])
end
set(gca, 'FontSize', textsize)

titlestr = append(sprintf('Total GF: %.2f MeV - %.2f MeV ', min_incident_energy, max_incident_energy), ...
    addin, sprintf(' %.0f Bins', bins));
%title(titlestr, 'FontSize', titlesize)
%title('Proton Total Geometric Factor', 'FontSize', titlesize-2);
%ylabel('Instrument Efficiency', 'FontSize', textsize)
ylabel('Geometric Factor (cm^2 sr)', 'FontSize', textsize)
xlabel('Incident Energy (MeV)', 'FontSize', textsize)

% Changes legend depending on the Sim_Type
%{
if sim_type == 0
    legend_entries = {'Theoretical Min', 'Theoretical Max', 'GEANT4 Cap'};
elseif sim_type == 1
    legend_entries = {'Theoretical Min', 'Theoretical Max', 'GEANT4 Sphere'};
end

legend(legend_entries, 'FontSize', titlesize, 'Location', 'southeast');
%}
hold off

% Saving the figure as a jpg then returning to the main directory
effsave = append('Total GF_', string(datetime("today")), '_', addin, '.png');
saveas(f1, effsave)

%% Geometric Factor by Energy Channel
geo_EC(geo_EC==0) = 1e-31;
back_bool = exist('geo_back_EC_combined','var');

% 1. Vectorized String Generation (Eliminates loops)
EngLegend_low  = string(compose('%.3f', energy_channels(:, 1)))';
EngLegend_high = string(compose('%.3f', energy_channels(:, 2)))';
EngLegend_high(end) = sprintf('%.0f', energy_channels(end, 2));

% 2. Create Plot Colors
Effplotcolor = plasma(size(geo_EC, 1)); 
if back_bool
    full_winter = winter(size(geo_EC, 1));
    Effplotcolor_back = full_winter((end - size(geo_back_EC_combined, 1) + 1):end, :);
end

channel_select = [];
channel_select = [5,10,15,20,25];
% 3. Default Channel Selection & Setup Figure
if isempty(channel_select), channel_select = 1:size(geo_EC, 1); end

f3 = figure('Position', [0 0 2800 1080]); hold on;

% 4. Consolidated Data Plotting
for i = 1:size(geo_EC, 1)
    is_sel = ismember(i, channel_select);
    c_prim = [0.7 0.7 0.7 0.7]; if is_sel, c_prim = Effplotcolor(i,:); end
    
    plot(energy_midpoints, geo_EC(i,:), 'Color', c_prim, 'LineWidth', 2);
    
    if back_bool
        if i <= 8, b_idx = 1; else, b_idx = i - 7; end
        c_back = [0.7 0.7 0.7 0.7]; if is_sel, c_back = Effplotcolor_back(b_idx,:); end
        plot(energy_midpoints, geo_back_EC_combined(b_idx,:), 'Color', c_back, 'LineWidth', 2);
    end
end

% 5. Apply Formatting (Condensed)
set(gca, 'FontSize', textsize, 'XScale', 'log', 'YScale', 'log');
ylabel('Geometric Factor (cm^2 sr)', 'FontSize', textsize);
xlabel('Incident Energy (MeV)', 'FontSize', textsize);
grid on;

if parttype == 0
    set(gca, 'XTick', [0,0.5,1,1.5,2,3,4,5,6,7,8,9,10]);
    xlim([0.5, max_incident_energy]); ylim([10^-4, 10^0]);
    xline([0.62, 5.71, 15.91, 30.89], '--', {'Be Penetration','D9 Penetration', 'W Penetration', 'W & Stack Penetration'}, ...
        'LineWidth', 1.5, 'FontSize', 16, 'LabelOrientation', 'aligned');
    title('Electron Energy Channels')
elseif parttype == 1
    set(gca, 'XTick', [0,10,20,30,40,50,60,80,100,120,140,160,200,300,400,500,600,800,1000]);
    xlim([10, max_incident_energy]); ylim([10^-4, 10^1.2]);
    xline([14.15, 52.29, 69.09, 99.64], '--', {'Be Penetration','D9 Triggering', 'W Penetration', 'W & Stack Penetration'}, ...
        'LineWidth', 1.5, 'FontSize', 16, 'LabelOrientation', 'aligned');
    title('Proton Energy Channels')
end
hold off;

% --- BEGIN CUSTOM RIGHT-SIDE DUAL LEGEND ---
% 1. Securely grab the main plot and force it to shrink to 70% width
ax_main = gca;
set(ax_main, 'Position', [0.05, 0.11, 0.75, 0.81]); 

% 2. Create the legend axes in the newly opened space on the right
% Starts at 0.78, leaving an 5% physical gap between the plot and the legend
ax_leg = axes('Position', [0.76, 0.11, 0.21, 0.81], 'Visible', 'off');
hold(ax_leg, 'on'); xlim(ax_leg, [0, 1]); ylim(ax_leg, [0, 1]);

N_leg = length(channel_select);
dy = min(0.03, 0.91 / max(1, N_leg - 1)); 
y_rows = 0.93 - (0:(N_leg-1)) * dy;

% Adjusted horizontal spacing to fit the newly sized legend box
x_text = 0.50;        % The absolute center of the legend
box_offset = 0.22;

x_range = x_text - box_offset; 
x_pen   = x_text + box_offset;

% Draw Headers
text(x_range, 0.98, 'Range', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', textsize-10);
text(x_text, 0.98, 'Energy (MeV)', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', textsize-10);
if back_bool
    text(x_pen, 0.98, 'Penetrating', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', textsize-10); 
end

y_combined = []; 
for i = 1:N_leg
    idx = channel_select(i); 
    
    if idx <= size(geo_EC, 1)
        plot(x_range, y_rows(i), 'sk', 'MarkerSize', 14, 'MarkerFaceColor', Effplotcolor(idx,:));
        
        fs = textsize - 4;
        tl = text(x_text - 0.04, y_rows(i), EngLegend_low(idx), 'HorizontalAlignment', 'right', 'FontSize', fs);
        tm = text(x_text, y_rows(i), '-', 'HorizontalAlignment', 'center', 'FontSize', fs);
        tr = text(x_text + 0.04, y_rows(i), EngLegend_high(idx), 'HorizontalAlignment', 'left', 'FontSize', fs);
            
        % Dynamic font shrinker to protect boundaries
        while (tl.Extent(1) < (x_range + 0.08) || (tr.Extent(1) + tr.Extent(3)) > (x_pen - 0.08)) && fs > 6
            fs = fs - 0.5; set([tl, tm, tr], 'FontSize', fs); drawnow limitrate;
        end
    end
    
    if back_bool
        if idx <= 8
            y_combined(end+1) = y_rows(i);
        else
            plot(x_pen, y_rows(i), 'sk', 'MarkerSize', 14, 'MarkerFaceColor', Effplotcolor_back(idx-7,:));
        end
    end
end

if back_bool && ~isempty(y_combined)
    plot(x_pen, mean(y_combined), 'sk', 'MarkerSize', 14, 'MarkerFaceColor', Effplotcolor_back(1,:));
end
% --- END CUSTOM RIGHT-SIDE DUAL LEGEND ---

saveas(f3, append('Combined_Geometric_Factor_Whole_by_EC_', string(datetime("today")), addin, '.png'));

%% Geometric Factor by Energy Channel (range particles)
geo_EC(geo_EC==0) = 1e-31;

% 1. Vectorized String Generation (Eliminates loops)
EngLegend_low  = compose('%.3f', energy_channels(:, 1))';
EngLegend_high = compose('%.3f', energy_channels(:, 2))';

% 2. Create Plot Colors
Effplotcolor = plasma(size(geo_EC, 1)); 

% 3. Default Channel Selection & Setup Figure
if ~exist('channel_select','var'), channel_select = 1:size(geo_EC, 1); end

f4 = figure('Position', [0 0 2800 1080]); hold on;

% 4. Consolidated Data Plotting
for i = 1:size(geo_EC, 1)
    is_sel = ismember(i, channel_select);
    c_prim = [0.7 0.7 0.7 0.7]; if is_sel, c_prim = Effplotcolor(i,:); end
    
    plot(energy_midpoints, geo_EC(i,:), 'Color', c_prim, 'LineWidth', 2);
end

% 5. Apply Formatting (Condensed)
set(gca, 'FontSize', textsize, 'XScale', 'log', 'YScale', 'log');
ylabel('Geometric Factor (cm^2 sr)', 'FontSize', textsize);
xlabel('Incident Energy (MeV)', 'FontSize', textsize);
grid on;

if parttype == 0
    set(gca, 'XTick', [0,0.5,1,1.5,2,3,4,5,6,7,8,9,10]);
    xlim([0.5, 10]); ylim([10^-4, 10^0]);
    xline([0.62, 5.71, 15.91, 30.89], '--', {'Be Penetration','D9 Penetration', 'W Penetration', 'W & Stack Penetration'}, ...
        'LineWidth', 1.5, 'FontSize', 16, 'LabelOrientation', 'aligned');
elseif parttype == 1
    set(gca, 'XTick', [0,10,20,30,40,50,60,80,100,120,140,160,200,300,400,500,600,800,1000]);
    xlim([10, max_incident_energy]); ylim([10^-4, 10^1.2]);
    xline([14.15, 52.29, 69.09, 99.64], '--', {'Be Penetration','D9 Triggering', 'W Penetration', 'W & Stack Penetration'}, ...
        'LineWidth', 1.5, 'FontSize', 16, 'LabelOrientation', 'aligned');
end
hold off;

% --- BEGIN CUSTOM RIGHT-SIDE SINGLE LEGEND ---
% 1. Securely grab the main plot and force it to shrink
ax_main = gca;
set(ax_main, 'Position', [0.05, 0.11, 0.75, 0.81]); 

% 2. Create the legend axes
ax_leg = axes('Position', [0.8, 0.11, 0.21, 0.81], 'Visible', 'off');
hold(ax_leg, 'on'); xlim(ax_leg, [0, 1]); ylim(ax_leg, [0, 1]);

N_leg = length(channel_select);
dy = min(0.03, 0.91 / max(1, N_leg - 1)); 
y_rows = 0.93 - (0:(N_leg-1)) * dy;

% Adjusted horizontal spacing for a 2-column layout
x_text = 0.5;
box_offset = 0.30;
x_rng = x_text - box_offset; 

% Draw Headers
text(x_rng, 0.98, 'Range', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', textsize-6);
text(x_text, 0.98, 'Energy (MeV)', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', textsize-6);

for i = 1:N_leg
    idx = channel_select(i); 
    
    plot(x_rng, y_rows(i), 'sk', 'MarkerSize', 14, 'MarkerFaceColor', Effplotcolor(idx,:));
    
    fs = textsize - 12;
    tl = text(x_text - 0.04, y_rows(i), EngLegend_low(idx), 'HorizontalAlignment', 'right', 'FontSize', fs);
    tm = text(x_text, y_rows(i), '-', 'HorizontalAlignment', 'center', 'FontSize', fs);
    tr = text(x_text + 0.04, y_rows(i), EngLegend_high(idx), 'HorizontalAlignment', 'left', 'FontSize', fs);
        
    % Dynamic font shrinker to protect boundaries
    left_limit = x_rng + 0.07;
    right_limit = 0.95; % No right box, so limit is just the edge of the legend space
    
    while (tl.Extent(1) < left_limit || (tr.Extent(1) + tr.Extent(3)) > right_limit) && fs > 6
        fs = fs - 0.5; set([tl, tm, tr], 'FontSize', fs); drawnow limitrate;
    end
end
% --- END CUSTOM RIGHT-SIDE SINGLE LEGEND ---

% Save the figure
effsave_rng = append('Range_Geometric_Factor_Whole_by_EC_', string(datetime("today")), addin, '.png');
saveas(f4, effsave_rng);

%% Geometric Factor by Energy Channel (penetrating particles)
if parttype == 1
    geo_back_EC_combined(geo_back_EC_combined==0) = 1e-31;
    
    % 1. Vectorized String Generation (Eliminates loops)
    % Uses the COMBINED energy channels matrix directly
    EngLegend_back_low  = compose('%.3f', energy_channels_combined(:, 1))';
    EngLegend_back_high = compose('%.3f', energy_channels_combined(:, 2))';
    
    % 2. Create Plot Colors (Shifted to match the combined channels)
    full_winter = winter(size(energy_channels, 1)); 
    Effplotcolor_back = full_winter((end - size(geo_back_EC_combined, 1) + 1):end, :);
    
    % 3. Map Original Channel Selections to the Combined Matrix Indices
    if exist('channel_select','var')
        back_select = [];
        for c = channel_select
            if c <= 8
                back_select(end+1) = 1;
            else
                back_select(end+1) = c - 7;
            end
        end
        back_select = unique(back_select);
    else
        back_select = 1:size(geo_back_EC_combined, 1);
    end
    
    % Setup Figure
    f5 = figure('Position', [0 0 2800 1080]); hold on;
    
    % 4. Consolidated Data Plotting using 1-to-1 mapped indices
    for back_idx = 1:size(geo_back_EC_combined, 1)
        is_sel = ismember(back_idx, back_select);
        c_back = [0.7 0.7 0.7 0.7]; if is_sel, c_back = Effplotcolor_back(back_idx, :); end
        
        plot(energy_midpoints, geo_back_EC_combined(back_idx,:), 'Color', c_back, 'LineWidth', 2);
    end
    
    % 5. Apply Formatting (Condensed)
    set(gca, 'FontSize', textsize, 'XScale', 'log', 'YScale', 'log');
    ylabel('Geometric Factor (cm^2 sr)', 'FontSize', textsize);
    xlabel('Incident Energy (MeV)', 'FontSize', textsize);
    grid on;
    
    if parttype == 0
        set(gca, 'XTick', [0,0.5,1,1.5,2,3,4,5,6,7,8,9,10]);
        xlim([0.5, max_incident_energy]); ylim([10^-4, 10^0]);
        xline([0.62, 5.71, 15.91, 30.89], '--', {'Be Penetration','D9 Penetration', 'W Penetration', 'W & Stack Penetration'}, ...
            'LineWidth', 1.5, 'FontSize', 16, 'LabelOrientation', 'aligned');
    elseif parttype == 1
        set(gca, 'XTick', [0,10,20,30,40,50,60,80,100,120,140,160,200,300,400,500,600,800,1000]);
        xlim([10, max_incident_energy]); ylim([10^-4, 10^1.2]);
        xline([14.15, 52.29, 69.09, 99.64], '--', {'Be Penetration','D9 Triggering', 'W Penetration', 'W & Stack Penetration'}, ...
            'LineWidth', 1.5, 'FontSize', 16, 'LabelOrientation', 'aligned');
    end
    hold off;
    
    % --- BEGIN CUSTOM RIGHT-SIDE SINGLE LEGEND ---
    % 1. Securely grab the main plot and force it to shrink
    ax_main = gca;
    set(ax_main, 'Position', [0.05, 0.11, 0.70, 0.81]); 
    
    % 2. Create the legend axes
    ax_leg = axes('Position', [0.78, 0.11, 0.21, 0.81], 'Visible', 'off');
    hold(ax_leg, 'on'); xlim(ax_leg, [0, 1]); ylim(ax_leg, [0, 1]);
    
    N_leg = length(back_select);
    dy = min(0.03, 0.91 / max(1, N_leg - 1)); 
    y_rows = 0.93 - (0:(N_leg-1)) * dy;
    
    % Adjusted horizontal spacing for a 2-column layout
    x_text = 0.35;        % Shifted left to balance the missing Range column
    box_offset = 0.30;    % Distance between the Energy text and Penetrating box
    x_pen = x_text + box_offset; 
    
    % Draw Headers
    text(x_text, 0.98, 'Energy (MeV)', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', textsize-6);
    text(x_pen, 0.98, 'Penetrating', 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', textsize-6);
    
    for i = 1:N_leg
        idx = back_select(i); 
        
        plot(x_pen, y_rows(i), 'sk', 'MarkerSize', 14, 'MarkerFaceColor', Effplotcolor_back(idx,:));
        
        fs = textsize - 12;
        tl = text(x_text - 0.04, y_rows(i), EngLegend_back_low(idx), 'HorizontalAlignment', 'right', 'FontSize', fs);
        tm = text(x_text, y_rows(i), '-', 'HorizontalAlignment', 'center', 'FontSize', fs);
        tr = text(x_text + 0.04, y_rows(i), EngLegend_back_high(idx), 'HorizontalAlignment', 'left', 'FontSize', fs);
            
        % Dynamic font shrinker to protect boundaries
        left_limit = 0.05; % Protect text from falling off the left edge of the legend space
        right_limit = x_pen - 0.07;
        
        while (tl.Extent(1) < left_limit || (tr.Extent(1) + tr.Extent(3)) > right_limit) && fs > 6
            fs = fs - 0.5; set([tl, tm, tr], 'FontSize', fs); drawnow limitrate;
        end
    end
    % --- END CUSTOM RIGHT-SIDE SINGLE LEGEND ---
    
    % Save the figure
    effsave_back = append('Penetrating_Geometric_Factor_Whole_by_EC_', string(datetime("today")), addin, '.png');
    saveas(f5, effsave_back);
end

%% Save Total Geometric Factor and Geofactor by EC
%
% Save total geometric factor for range particles
filename = sprintf('electron_FS_31000_v2_range_totalGF.txt');
fileID = fopen(filename, 'w');
for i = 1:length(energy_midpoints)
fprintf(fileID,'%.6f %.6f\n',energy_midpoints(i),total_geo(i));
end
fclose(fileID);

% Save total geometric factor for penetrating particles
filename = sprintf('electron_FS_31000_v2_pen_totalGF.txt');
fileID = fopen(filename, 'w');
for i = 1:length(energy_midpoints)
fprintf(fileID,'%.6f %.6f\n',energy_midpoints(i),back_total_geo(i));
end
fclose(fileID);

% Save geometric factor by energy channel for range particles
filename = sprintf('electron_FS_31000_v2_range_GFbyEC.txt');
fileID = fopen(filename, 'w');
% 1. Determine the variable number of rows in geo_EC
num_geo_rows = size(geo_EC, 1);
% 2. Dynamically build the format string
% This creates one '%.6f' for the energy, plus ' %.6f' for every row in geo_EC
formatSpec = ['%.6f', repmat(' %.6f', 1, num_geo_rows), '\n'];
% 3. Execute the print loop
for i = 1:length(energy_midpoints)
    % geo_EC(:,i)' is transposed to ensure it feeds into fprintf as a row vector
    fprintf(fileID, formatSpec, energy_midpoints(i), geo_EC(:,i)');
end
fclose(fileID);

% Save geometric factor by energy channel for penetrating particles
filename = sprintf('electron_FS_31000_v2_pen_GFbyEC.txt');
fileID = fopen(filename, 'w');
% 1. Determine the variable number of rows in geo_EC
num_geo_rows = size(geo_back_EC, 1);
% 2. Dynamically build the format string
% This creates one '%.6f' for the energy, plus ' %.6f' for every row in geo_EC
formatSpec = ['%.6f', repmat(' %.6f', 1, num_geo_rows), '\n'];
% 3. Execute the print loop
for i = 1:length(energy_midpoints)
    % geo_EC(:,i)' is transposed to ensure it feeds into fprintf as a row vector
    fprintf(fileID, formatSpec, energy_midpoints(i), geo_back_EC(:,i)');
end
fclose(fileID);
%


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
                    M_energy_bin = 2 .* M_energy_bin / (1 - cosd(source_angle));
                    
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
                    geo_EC(channel,:) = hits_log(channel,:) ./ M_energy_bin * (4 * (pi^2) * (r_source^2));
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
cd 'D:\HERT_Drive\MATLAB Main\Result'