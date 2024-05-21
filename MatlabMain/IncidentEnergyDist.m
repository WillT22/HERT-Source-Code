% Resets all variables and values in MATLAB
clear all;
close all;
clc;
addpath 'E:\HERT_Drive\MATLAB Main'

files = dir('E:\HERT_Drive\Matlab Main\Result\Proton\RAW_data');
% Grabs all the names of the files in a vector (nx1 matrix)
filenames = {files.name};
filenames = filenames(3:end);
file_name = filenames{1};
beam_number = str2double(regexp(file_name, '\d+', 'match', 'once'));

Einc = zeros(beam_number,size(filenames,2));

for f = 1:size(filenames,2)
    % Opens the file to be processed
    fide = fopen(filenames{f}, 'r');
    % Sorts each line in the .txt file into a cell array
    header = fgetl(fide);  % Read 11 strings directly into a cell array
    data = textscan(fide, '%f %f %f %f %f %f %f %f %f %f', 'Delimiter','');  % Skip header
    % Closes the file
    fclose(fide);
                    
    % Separates data into incident energy and detected energy
    Einc(:,f) = data{1};
end

Einc = reshape(Einc,1,size(Einc,1)*size(Einc,2));
Einc_sort = sort(Einc);

binWidth = 0.1;  % Desired bin width
lowerBound = 10;   % Lower bound of data
upperBound = 80;  % Upper bound of data

binEdges = lowerBound:binWidth:upperBound;  % Create bin edges

% Create the histogram
histogramData = histogram(Einc_sort, 'BinEdges', binEdges);

% Normalize the frequencies to percentages
totalParticles = length(Einc_sort);  % Assuming Einc_sort contains all particles
normalizedFrequencies = histogramData.Values / totalParticles * 100;

% Plot the normalized histogram
bar(binEdges(1:end-1), normalizedFrequencies);

% Customize the plot
title('Histogram with Percentage of Total Particles');
xlabel('Data Value');
ylabel('Percentage of Total Particles (%)');
