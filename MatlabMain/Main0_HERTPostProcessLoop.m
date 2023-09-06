%% Skyler Krantz
% Edited by Will Teague 
% Updated: Sept. 6th, 2023
% Post Processing for GEANT4

%Resets all variables and values in MATLAB
clear all;
close all;
clc;
addpath 'E:\HERT\Matlab Main\'
%% Reads in all files names

cd 'E:\HERT\Matlab Main\Result\'; %Changes address to where results are stored
%read files from ./Result folder stores into 1*C array
%inputfolder = 'E:\HERT\Matlab Main Share\Result';
outputfolder = 'E:\HERT\Matlab Main\Result\PostProcess';%location for where post process files should be saved

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

FolderChoice = menu('HERT Post Process Loop: Choose an Folder',subFolderNames{:});
inputfolder = subFolderNames{FolderChoice};
inputfolder = append(pwd,'\',inputfolder);
%Creates search string for result .txt files
inputfiles = append(inputfolder,'/*.txt');

%lists all .txt files in the input folder
list = dir(inputfiles);

%Grabs all the names of the files in vector (nx1 matrix)
list_fileNames = {list.name};

%Gets number of columns in file names
%indicating number of files that can be loaded
C = size(list_fileNames,2);

choice = menu('HERT Post Process Loop: Choose an option', 'Exit Program', 'Load one file','Load all files');

while choice~=1
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
                    fprintf('Number of files loaded: %.0f \n',length(list_fileNames))
            end
    
         case 4
            %One file selected
            if size(filename,1)==1
                disp('Run HERTPostProcessWhole.m Once');
                
                Output = HERTPostProcessWhole(filename,inputfolder,outputfolder);
            % More than one file selected
            elseif size(filename,1) > 1
                disp('Start to loop HERTPostProcessWhole.m');

                %Nested For loops to create final matrix 1 and 2
                for i = 1: C 
                    tic

                    Output = HERTPostProcessWhole(filename{i},inputfolder,outputfolder);

                    toc
                end
    
            end
    end
    choice = menu('HERT Post Process Loop:Choose an option', 'Exit Program', 'Load one file','Load all files','Start Run');
end