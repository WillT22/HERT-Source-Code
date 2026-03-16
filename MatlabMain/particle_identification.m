% Requires data from Main1
min_dep_energy = 35;
problem_particles = particles(M_hit_dep>M_energy_beam-0.01 & M_energy_beam<40,:);

% Choosing select / random runs only
selected_runs = [];
%selected_runs = [2540,4614];
%random_runs = 10;
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
cd 'D:\HERT_Drive\Matlab Main\Result\Proton_FS\raw_data\'

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
fileID = fopen('proton_fullpen.csv','w');
for r = 1:size(problem_particles,1)
    fprintf(fileID,['%10.6f, %10.0f, %10.0f, %10.6f, %10.6f, %10.6f, %10.6f,' ...
        ' %10.6f, %10.6f, %10.6f, %10.6f, %10.6f,\n'],problem_particles(r,:));
end
fclose(fileID);
%

%% Plot Edep v Einc
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
