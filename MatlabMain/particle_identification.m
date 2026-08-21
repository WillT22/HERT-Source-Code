% Requires data from Main1
min_dep_energy = 0;
max_dep_energy = 5.5;
min_energy_beam = 35;
max_energy_beam = 38;
problem_particles = particles(M_hit_dep>min_dep_energy & M_hit_dep<max_dep_energy...
                                & M_energy_beam>min_energy_beam &M_energy_beam<max_energy_beam,:);

% Choosing select / random runs only
selected_runs = [];
%selected_runs = [2540,4614];
random_runs = [];
random_runs = 20;
if ~isempty(selected_runs)
    problem_particles = problem_particles(ismember(problem_particles(:,2), selected_runs),:);
elseif isempty(selected_runs) && ~isempty(random_runs)
    unique_runs = unique(problem_particles(:,2));
    unique_runs_limit = unique_runs(unique_runs>2000 & unique_runs<10000);
    rand_idx = randperm(length(unique_runs_limit),random_runs);
    selected_runs = unique_runs_limit(rand_idx);
    selected_runs = sort(selected_runs);
    problem_particles = problem_particles(ismember(problem_particles(:,2), selected_runs),:);
end
%
%cd 'D:\HERT_Drive\MATLAB Main\Result'; % Main Result Directory
cd '\\arc.auburn.edu\lab\hert\Will Results\Proton_FS\raw_data'

for i = 1:length(problem_particles(:,2))
    run_number = problem_particles(i,2);
    filename = sprintf('HERT_CADoutput_proton_1000000_Run%d.txt', run_number);
    
    if exist(filename, 'file') == 2
        fid = fopen(filename, 'r');
        fprintf('Processing Run %d\n', run_number);
        if fid == -1
            warning('Could not open file: %s', filename);
            continue;
        end
        
        % Read and process file here
        data = textscan(fid, '%f %f %f %f %f %f %f %f %f %f', 'Delimiter','','HeaderLines',1);  % Skip header
        
        fclose(fid);

        Einc = data{1};
        Detector_Energy = cell2mat(data(2:end));

        problem_tempE = problem_particles(problem_particles(:,2) == run_number, 1);  % Get all energies for this run
        for j = 1:length(problem_tempE)
            current_energy = problem_tempE(j);
            particle_numbers = find(Einc == current_energy);
            valid_particles = particle_numbers(sum(Detector_Energy(particle_numbers,:),2) > min_dep_energy);
            
            if ~isempty(valid_particles)
                matching_row = find(problem_particles(:,1) == current_energy & problem_particles(:,2) == run_number);
                problem_particles(matching_row, 3) = valid_particles(1:length(matching_row));  % Take first valid particle if multiple
                problem_particles(matching_row, 4:12) = Detector_Energy(valid_particles(1:length(matching_row)),:);
            end
        end

    else
        warning('File does not exist: %s', filename);
    end
end

%
cd ../../
clear fileID
fileID = fopen('proton_colpen.csv','w');
for r = 1:size(problem_particles,1)
    fprintf(fileID,['%10.6f, %10.0f, %10.0f, %10.6f, %10.6f, %10.6f, %10.6f,' ...
        ' %10.6f, %10.6f, %10.6f, %10.6f, %10.6f,\n'],problem_particles(r,:));
end
fclose(fileID);
%

%% Plot Edep v Einc
%{
% Prepare data for plotting
M_energy_beam_reduced = M_energy_beam(M_energy_beam < 100);
M_hit_dep_reduced = M_hit_dep(M_energy_beam < 100);
random_indices = randperm(length(M_energy_beam_reduced), 500000);
M_energy_beam_subset = M_energy_beam_reduced(random_indices);
M_hit_dep_subset = M_hit_dep_reduced(random_indices);

figure; % Open a new figure window

% X-axis: M_energy_beam (Incident Beam Energy)
% Y-axis: M_hit_dep (Measured Deposited Energy)
scatter(M_energy_beam_subset, M_hit_dep_subset, 10, 'filled', 'MarkerFaceColor', [0 0.4470 0.7410]); % Marker size reduced to 1 for large N

% --- 5. Enhance Plot Aesthetics ---
hold on;
% Calculate max_val based on the full range of the data (not just the subset)
max_val = max([M_energy_beam, M_hit_dep]) * 1.1; 

% Add a 1:1 line for visual comparison
plot([0.1, max_val], [0.1, max_val], 'k--', 'LineWidth', 1.5, 'DisplayName', '1:1 Line'); 
hold off;

title('Deposited Energy vs. Incident Beam Energy', 'FontSize', 18);
xlabel('Incident Beam Energy (MeV)', 'FontSize', 16); 
ylabel('Measured Deposited Energy (MeV)', 'FontSize', 16); 
legend('Random Subset Data', '1:1 Line', 'Location', 'northwest','FontSize', 16);

grid on; 
xlim([10, max(M_energy_beam_subset) * 1.1]);
ylim([1, max(M_hit_dep) * 1.1]);
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'XTick', [0,10,20,30,40,50,60,80,100,130,160,200,300,400]);
set(gca, 'YTick', [0,10,20,30,40,50,60,70,80]);
set(gca, 'FontSize', 16)

disp('');
%}

%%
%{
% --- Target Parameters ---
target_Einc = 18.8241;
min_dep_energy = 0.1;
tolerance = 1e-4; % Crucial for floating-point comparison

% Navigate to raw data directory
% cd 'D:\HERT_Drive\MATLAB Main\Result'; % Main Result Directory
cd '\\arc.auburn.edu\lab\hert\Will Results\Proton_FS\raw_data'

% Get list of all matching text files in the directory
file_list = dir('HERT_CADoutput_proton_1000000_Run9*.txt');
fprintf('Combing through %d files for E_inc = %.4f...\n', length(file_list), target_Einc);

% Initialize array to store results: 
% [Run Number, Particle Row, Total Deposited Energy, D1_energy, D2_energy...]
found_cases = []; 

total_files = length(file_list);
for i = 1:total_files
    filename = file_list(i).name;
    fprintf('Scanning file %d of %d: %s\n', i, total_files, filename);
    
    % Extract run number directly from the filename string
    run_number = sscanf(filename, 'HERT_CADoutput_proton_1000000_Run%d.txt');
    
    fid = fopen(filename, 'r');
    if fid == -1
        warning('Could not open file: %s', filename);
        continue;
    end
    
    % Read and process file
    data = textscan(fid, '%f %f %f %f %f %f %f %f %f %f', 'Delimiter', '', 'HeaderLines', 1);
    fclose(fid);
    
    Einc = data{1};
    Detector_Energy = cell2mat(data(2:end));
    total_dep = sum(Detector_Energy, 2);
    
    % Find rows matching the target incident energy AND minimum deposited energy
    % abs() is used to prevent floating-point mismatch errors
    matching_rows = find(abs(Einc - target_Einc) < tolerance & total_dep > min_dep_energy);
    
    if ~isempty(matching_rows)
        for m = 1:length(matching_rows)
            idx = matching_rows(m);
            fprintf('>>> MATCH FOUND! Run: %d | Row: %d | Total Dep: %.4f MeV\n', ...
                run_number, idx, total_dep(idx));
            
            % Store the result for extraction/printing
            found_cases = [found_cases; run_number, idx, total_dep(idx), Detector_Energy(idx, :)];
        end
    end
end

cd ../../

% Output results to a specific target file
if ~isempty(found_cases)
    fileID = fopen('target_particle_18_8241.csv','w');
    fprintf(fileID, 'Run_Num, Row_Idx, Total_Dep, D1, D2, D3, D4, D5, D6, D7, D8, D9\n');
    for r = 1:size(found_cases, 1)
        fprintf(fileID, '%d, %d, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f\n', found_cases(r,:));
    end
    fclose(fileID);
    fprintf('Target particles saved to target_particle_18_8241.csv\n');
else
    fprintf('No particles matching E_inc = %.4f with Dep > %d MeV were found.\n', target_Einc, min_dep_energy);
end
%}