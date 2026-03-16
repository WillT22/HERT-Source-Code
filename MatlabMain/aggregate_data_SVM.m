% Compiling all hits into one file 
clear all;
close all;
clc;
addpath 'D:\HERT_Drive\MATLAB Main'

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

% Creates a search string for result .txt files
file_prefix = 'PostProcessHERT';
file_suffix = '.txt';
inputfiles = append(inputfolder, '/', file_prefix, '*', file_suffix);
% Lists all .txt files in the Result folder
list = dir(inputfiles);
% Grabs all the names of the files in a vector (nx1 matrix)
list_fileNames = {list.name};

num_files_to_select = 200;
total_available_files = length(list_fileNames);

if total_available_files > num_files_to_select
    % Generate 1000 unique random integers between 1 and total_available_files
    random_indices = randperm(total_available_files, num_files_to_select);
    % Subset the file names list using the random indices
    list_fileNames = list_fileNames(random_indices);
    fprintf('Randomly selected %d files out of %d available.\n', num_files_to_select, total_available_files);
else
    fprintf('Warning: Only %d files found. Processing all available files instead of %d.\n', total_available_files, num_files_to_select);
end

% Gets the number of rows and columns in file names. Columns will indicate
% the number of files that can be loaded
file_number = length(list_fileNames);

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
    run_number = str2double(regexp(file_name{i}, '\d+_Run(\d+)', 'tokens', 'once'));
    fprintf('Starting Run %d \n', run_number)
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
aggregate = round(aggregate,6);

tot_Edep = sum(NumEnergyDeposit);
fprintf('Total Energy Depositing Particles %g \n',tot_Edep);

cd 'Aggregate Data'\
fid = fopen('Aggregate Proton_FS Data.txt', 'wt');
fprintf(fid, '%s \n', header);
fprintf(fid, '%9.6g           %9.6g %9.6g %9.6g %9.6g %9.6g %9.6g %9.6g %9.6g %9.6g \n', aggregate');
fclose(fid);

% Return to MATLAB Main
cd ..\..