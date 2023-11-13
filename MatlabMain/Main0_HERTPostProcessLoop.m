%% Authored by Skyler Krantz
% Edited by Will Teague 
% Updated: Nov. 13th, 2023
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

% Gets number of columns in file names indicating the number of files that can be loaded
C = size(list_fileNames, 2);

% User chooses an option from the menu
choice = menu('HERT Post Process Loop: Choose an option', 'Exit Program', 'Load one file', 'Load all files');

while choice ~= 1
    switch choice
        case 0
            disp('Error - please choose one of the options.')

        % Load One File
        case 2
            % Displays all the files that can be loaded
            txt_file_choice = menu('Choose a file', list_fileNames{:});

            % If menu button is closed, it will recycle to the initial menu
            if txt_file_choice == 0
                disp('Please select a file')
            % If the user selects a file, it will change filename to reference that one file.
            elseif txt_file_choice > 0
                for i = 1:C
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
            filename = list_fileNames';
            fprintf('Number of files loaded: %.0f \n', length(list_fileNames))

         case 4
            % One file selected
            if size(filename, 1) == 1
                disp('Run HERTPostProcessWhole.m Once');
                
                % Calls HERTPostProcessWhole function for one file
                Output = HERTPostProcessWhole(filename, inputfolder, outputfolder);
            
            % More than one file selected
            elseif size(filename, 1) > 1
                disp('Start to loop HERTPostProcessWhole.m');

                % Iterating over each file
                for i = 1:C 
                    tic

                    % Calls HERTPostProcessWhole function in a loop for each file
                    Output = HERTPostProcessWhole(filename{i}, inputfolder, outputfolder);

                    toc
                end
            end
    end
    
    % User chooses another option from the menu
    choice = menu('HERT Post Process Loop: Choose an option', 'Exit Program', 'Load one file', 'Load all files', 'Start Run');
end
