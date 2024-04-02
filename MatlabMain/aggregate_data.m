% Compiling all hits into one file 
clear all;
close all;
clc;
addpath 'E:\HERT_Drive\MATLAB Main'

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

% Creates a search string for result .txt files
inputfiles = append(inputfolder, '/*.txt');
% Lists all .txt files in the Result folder
list = dir(inputfiles);
% Grabs all the names of the files in a vector (nx1 matrix)
list_fileNames = {list.name};
% Gets the number of rows and columns in file names. Columns will indicate
% the number of files that can be loaded
file_number = size(list_fileNames, 2);

% Load all files
file_name = list_fileNames;
fprintf('Number of files loaded: %.0f\n', length(list_fileNames));
% Change to the input folder
cd(inputfolder);

NumEnergyDeposit = ones(1,file_number);
M_hit_dep = [];
M_Einc = [];

for i = 1:file_number
% Read the data
    fide = fopen(file_name{i}, 'r');
    % Sorts each line in the .txt file into a cell array
    NumEnergyDeposit_cell = textscan(fide, 'Sims with Energy Deposited: %f', 'Delimiter','');
    header = fgetl(fide);
    deposit_data = textscan(fide, '%f %f %f %f %f %f %f %f %f %f', 'Delimiter','');  % Skip header
    % Closes the file
    fclose(fide);

    % Extract relevant information from the data
    NumEnergyDeposit(i) = NumEnergyDeposit_cell{1, 1};
    Einc = deposit_data{1};
    Detector_Energy = cell2mat(deposit_data(2:end));

    M_Einc = [M_Einc; Einc];
    M_hit_dep = [M_hit_dep; Detector_Energy];
end

aggregate = [M_Einc, M_hit_dep];

tot_Edep = sum(NumEnergyDeposit);
fprintf('Total Energy Depositing Particles %g \n',tot_Edep);

cd 'Aggregate Data'\
fid = fopen('Aggregate Electron Data 1-100.txt', 'wt');
fprintf(fid, '%s \n', header);
fprintf(fid, '%9.6g          %10.6g%10.6g%10.6g%10.6g%10.6g%10.6g%10.6g%10.6g%10.6g\n', aggregate');
fclose(fid);

% Return to MATLAB Main
cd ..\..