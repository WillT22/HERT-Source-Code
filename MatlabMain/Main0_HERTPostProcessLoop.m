%% Authored by Skyler Krantz
% Edited by Will Teague 
% Updated: Feb. 13th, 2024
% Post Processing Itteration Loop for GEANT4

% Resets all variables and values in MATLAB
clear all;
close all;
clc;
addpath 'D:\HERT_Drive\Matlab Main\'
%% Reads in all files names

% Changes address to where results are stored
cd 'D:\HERT_Drive\Matlab Main\Result\'; 

% Location for where post process files should be saved
outputfolder = 'D:\HERT_Drive\Matlab Main\Result\PostProcess'; 

textsize = 20;
titlesize = 28;

% Get Folder Names for User
topLevelFolder = pwd; % Current folder
files = dir(topLevelFolder);
dirFlags = [files.isdir];
subFolders = files(dirFlags);
subFolderNames = {subFolders(3:end).name}; % Extracts only folder names

% User chooses a folder from the list
FolderChoice = menu('HERT Post Process Loop: Choose a Folder', subFolderNames{:});
inputfolder = subFolderNames{FolderChoice};
inputfolder = append(pwd, '\', inputfolder);

% Creates search string for result .txt files
inputfiles = append(inputfolder, '/*.txt');

% Lists all .txt files in the input folder
list = dir(inputfiles);

% Grabs all the names of the files in vector (nx1 matrix)
list_fileNames = {list.name};
number_files = size(list_fileNames, 2); % Gets number of files that can be loaded

% User chooses an option from the menu
choice = menu('HERT Post Process Loop: Choose an option', 'Exit Program', 'Choose file', 'Load all files');

all_files_loaded = false;

while choice ~= 1 
    switch choice
        case 0
            disp('Error - please choose one of the options.')

        % Choose File
        case 2
            % Displays all the files that can be loaded
            txt_file_choice = menu('Choose a file', list_fileNames{:});

            % If menu button is closed, it will recycle to the initial menu
            if txt_file_choice == 0
                disp('Please select a file')
            % If the user selects a file, it will change filename to reference that file.
            elseif txt_file_choice > 0
                for i = 1:number_files
                    switch txt_file_choice
                        case i
                            filename = list_fileNames{i};
                    end
                end
                % Shows the file that was loaded
                disp(filename);
            end

        % Load all files
        case 3
            filename = list_fileNames;
            fprintf('Number of files loaded: %.0f \n', length(list_fileNames))

        % Start Run
        case 4
            min_Energy = 0;
            max_Energy = 1000;
            energy_edges = min_Energy:0.01:max_Energy;
            energy_average = (energy_edges(1:end-1)+energy_edges(2:end))/2;
            
            Inc_hist = zeros(size(energy_average));
            for file_index = 1:number_files
                tic % finds elapsed time of reading each data file
                % processes each file
                Einc_out = HERTPostProcessWhole(filename{file_index}, inputfolder,outputfolder);
                Inc_hist = Inc_hist + histcounts(Einc_out,energy_edges);
                toc
            end
            
            f = figure;
            f.Position = [0 0 1920 1080];
            hold on
            scatter(energy_average,Inc_hist,'.');
            hold off
            
            xlabel('Energy (MeV)', 'FontSize', textsize); % Label the x-axis
            ylabel('Count', 'FontSize', textsize);         % Label the y-axis
            title('Distribution of Incident Energy'); % (Optional) Add a title to the plot
            
            % Saving the figure as a png then returning to the main directory
            effsave = append('Incident Energy Counts', string(datetime("today")), '.png');
            saveas(f, effsave)


            % Report all files are loaded and move on
            fprintf('All Files Processed\n')
    end % switch end
    
    % User chooses another option from the menu
    choice = menu('HERT Post Process Loop: Choose an option', 'Exit Program', 'Load one file', 'Load all files', 'Start Run');
end % while not Exit Program
