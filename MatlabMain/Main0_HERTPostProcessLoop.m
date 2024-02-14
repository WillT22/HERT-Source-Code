%% Authored by Skyler Krantz
% Edited by Will Teague 
% Updated: Feb. 13th, 2024
% Post Processing Itteration Loop for GEANT4

% Resets all variables and values in MATLAB
clear all;
close all;
clc;
addpath 'E:\HERT_Drive\Matlab Main\'
%% Reads in all files names

% Changes address to where results are stored
cd 'E:\HERT_Drive\Matlab Main\Result\'; 

% Location for where post process files should be saved
outputfolder = 'E:\HERT_Drive\Matlab Main\Result\PostProcess'; 

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
            for file_index = 1:number_files
                tic % finds elapsed time of reading each data file
                % processes each file
                HERTPostProcessWhole(filename{file_index}, inputfolder,outputfolder);
                toc
            end
            % Report all files are loaded and move on
            fprintf('All Files Processed\n')
    end % switch end
    
    % User chooses another option from the menu
    choice = menu('HERT Post Process Loop: Choose an option', 'Exit Program', 'Load one file', 'Load all files', 'Start Run');
end % while not Exit Program
