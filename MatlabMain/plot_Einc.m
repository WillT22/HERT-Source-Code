% Resets all variables and values in MATLAB
clear;
close all;
clc;

% Add path to your functions
addpath('D:\HERT_Drive\Matlab Main\');

%% 1. Select Particle Type and Folders
% Select Particle Type Menu
parttype_choice = menu('Particle Type', 'Electron', 'Proton', 'Helium');
if parttype_choice == 0
    disp('Selection canceled. Exiting.');
    return;
end

% Map MATLAB's 1-based menu choice (1, 2, 3) to your logic (0, 1, 2)
parttype = parttype_choice - 1;

% Determine the specific string identifier for the file search
if parttype == 0
    particle_str = 'electron';
elseif parttype == 1
    particle_str = 'proton';
elseif parttype == 2
    % Change 'helium' to 'alpha' here if your GEANT4 output files use alpha instead
    particle_str = 'helium'; 
end

% Opens a standard UI folder dialog
disp('Waiting for folder selection...');
inputfolder = uigetdir(pwd, 'Select the folder containing the GEANT4 .txt result files');

if inputfolder == 0
    disp('Folder selection canceled. Exiting script.');
    return;
end

% Automatically create the PostProcess folder
outputfolder = fullfile(inputfolder, 'PostProcess_Histograms');
if ~exist(outputfolder, 'dir')
    mkdir(outputfolder);
end

% Find ONLY the specific text files matching the chosen particle type
search_pattern = sprintf('PostProcessHERT_CADoutput_%s_1000000_*.txt', particle_str);
inputfiles = fullfile(inputfolder, search_pattern);
list = dir(inputfiles);
list_fileNames = {list.name};
number_files = length(list_fileNames);

if number_files == 0
    fprintf('No matching %s .txt files found in %s\n', particle_str, inputfolder);
    return;
end

fprintf('Successfully found %d files to process.\n\n', number_files);

%% 2. Initialize Histogram Variables
textsize = 20;
titlesize = 28;

% Finds energy edges and midpoints based off of bins
if parttype == 0
    bins = 400;
    energy_min = 0.01;
    energy_max = 10;
    energy_edges = logspace(log10(energy_min),log10(energy_max),bins+1);
    energy_edges(end) = energy_max + 0.01;
elseif parttype == 1
    bins = 600;
    energy_min = 10;
    energy_max = 10000;
    energy_edges = logspace(log10(energy_min),log10(energy_max),bins+1);
elseif parttype == 2
    bins = 600;
    energy_min = 10;
    energy_max = 100000;
    energy_edges = logspace(log10(energy_min),log10(energy_max),bins+1);
end

% finding midpoints
energy_midpoints = (energy_edges(2:end) + energy_edges(1:end-1))/2;
bin_widths = diff(energy_edges);

% Pre-allocate the counts array (must be 1 element shorter than edges array)
Inc_hist = zeros(1, length(energy_edges) - 1);

% Initialize an empty list to store the names of corrupted files
bad_files = {};

%% 3. Process All Files
disp('Starting Batch Processing...');

for file_index = 1:number_files
    tic; % Start timer for this file
    current_file = list_fileNames{file_index};
    file_path = fullfile(inputfolder, current_file);
    
    % Use try...catch to prevent the script from crashing on bad/empty files
    try
        % 1. Open the file
        fid = fopen(file_path, 'r');
        
        % 2. Skip the top two header lines
        fgetl(fid); 
        fgetl(fid);
        
        % 3. Block 1: Read the incident energies for the "Hits"
        data_hits = textscan(fid, '%f %*[^\n]');
        
        % 4. Safely advance the file pointer past the middle text
        while ~feof(fid)
            line = fgetl(fid);
            if contains(line, 'Einc(MeV)')
                break; % Stop skipping once we hit the column header for the misses
            end
        end
        
        % 5. Block 2: Read the incident energies for the "Misses"
        data_misses = textscan(fid, '%f');
        
        % 6. Close the file immediately to free memory
        fclose(fid);
        
        % Combine both blocks into the final 1D array
        Einc_out = [data_hits{1}; data_misses{1}];
        
        % Bin the incident energies and add to the running total
        Inc_hist = Inc_hist + histcounts(Einc_out, energy_edges);
        
        % Print progress
        fprintf('Processed %d of %d: %s (%.2f s)\n', file_index, number_files, current_file, toc);
        
    catch ME
        % Safety check: If the script crashes mid-read, ensure the file isn't left locked open
        if exist('fid', 'var') && fid ~= -1
            fclose(fid);
        end
        % Log it and skip
        fprintf('*** ERROR SKIPPING FILE: %s ***\n', current_file);
        bad_files{end+1} = current_file; 
    end
end

%% 4. Error Reporting
if ~isempty(bad_files)
    fprintf('\n========================================\n');
    fprintf('SCAN COMPLETE: Found %d corrupted files.\n', length(bad_files));
    disp('The following files were skipped:');
    disp(bad_files');
    fprintf('========================================\n\n');
else
    fprintf('\nAll %d files processed successfully with no corruption.\n', number_files);
end

%% 5. Plot and Save Histogram (Publication Format)

% Setup Publication Typography & Physical Dimensions
textsize = 14;              
fig_width_in  = 7.0;        
fig_height_in = 2.8;        

f = figure('Units', 'inches', ...
           'Position', [1, 1, fig_width_in, fig_height_in], ...
           'PaperUnits', 'inches', ...
           'PaperPosition', [0, 0, fig_width_in, fig_height_in], ...
           'PaperPositionMode', 'auto'); 
ax = axes(f);
hold(ax, 'on');

% Calculate bin widths (MeV) and divide counts by the width to get Counts/MeV
bin_widths = diff(energy_edges);
Inc_hist_per_MeV = Inc_hist ./ bin_widths;

% Use the pre-binned histogram plotting method with the normalized counts
histogram(ax, 'BinEdges', energy_edges, 'BinCounts', Inc_hist_per_MeV, ...
    'FaceColor', 'b', 'EdgeColor', 'none');

% --- FORMATTING ---
set(ax, 'FontSize', textsize, 'Layer', 'top');
set(ax, 'XScale', 'log', 'YScale', 'log');

if parttype == 0
    set(gca, 'XTick', [0, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
    xlim([0.5, 10]); 
    xlabel('Kinetic Energy (MeV)', 'FontSize', textsize);

elseif parttype == 1
    set(ax, 'XTick', [10, 20, 30, 50, 100, 200, 300, 500, 1000, 2000, 5000, 10000]);
    xlim(ax, [10, 10000]);
    xlabel('Kinetic Energy (MeV)', 'FontSize', textsize);

elseif parttype == 2
    set(ax, 'XTick', [10, 20, 30, 50, 100, 200, 300, 500, 1000, 2000, 5000, 10000]);
    xlim(ax, [10, 10000]);
    xlabel('Kinetic Energy (MeV/nuc)', 'FontSize', textsize);
end

current_ylims = ylim(ax);
min_power = floor(log10(current_ylims(1)));
max_power = ceil(log10(current_ylims(2)));
set(ax, 'YTick', 10.^(min_power:max_power));
ylim(ax, [10.^(min_power),10.^(max_power)])

grid(ax, 'on'); 
ax.YMinorGrid = 'off';

ylabel(ax, 'Counts / MeV', 'FontSize', textsize);
xlabel(ax, 'Energy (MeV)', 'FontSize', textsize);
title(ax, 'Distribution of Incident Energy', 'FontSize', textsize);

hold(ax, 'off');

% --- SAVING MULTIPLE FORMATS ---
% Determine File Name
timestamp = string(datetime("today"));
base_filename = fullfile(outputfolder, sprintf('Incident_Energy_Counts_%s', timestamp));

% Save the interactive MATLAB figure file (.fig)
savefig(f, base_filename);
% Export High-Resolution PNG (300 DPI)
print(f, [base_filename, '.png'], '-dpng', '-r300');
fprintf('Publication figures successfully saved to:\n%s\n', outputfolder);