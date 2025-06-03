% Requires data from Main1
problem_particles = particles(particles(:,1) < 8, :);

% Choosing select / random runs only
%
selected_runs = [15001:15400];
%rand_idx = randperm(max(problem_particles(:,2)));
%selected_runs = rand_idx(1:5);
if ~isempty(selected_runs)
    problem_particles = problem_particles(ismember(problem_particles(:,2), selected_runs),:);
end
%
%cd 'E:\HERT_Drive\MATLAB Main\Result'; % Main Result Directory
cd 'A:\Will Results\Electron_FS\raw_data'

unique_runs = unique(problem_particles(:,2));
%unique_runs = unique(particles(:, 2));
problem_particles = zeros(length(unique_runs), size(particles, 2)); 

for i = 1:length(unique_runs)
    run_number = unique_runs(i);
    % to find the first hit particles in each file:
    run_particles = particles(particles(:, 2) == run_number, :); 
    if ~isempty(run_particles)
       problem_particles(i, 1:2) = run_particles(1, :); 
    end

    filename = sprintf('HERT_CADoutput_electron_1000000_Run%d.txt', run_number);
    
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
            valid_particles = particle_numbers(sum(Detector_Energy(particle_numbers,:),2) >= 0.1);
            
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
fileID = fopen('electron_first_final_rev2.csv','w');
for r = 1:size(problem_particles,1)
    fprintf(fileID,['%10.6f, %10.0f, %10.0f, %10.6f, %10.6f, %10.6f, %10.6f,' ...
        ' %10.6f, %10.6f, %10.6f, %10.6f, %10.6f,\n'],problem_particles(r,:));
end
fclose(fileID);
%