%% ERROR ANALYSIS FOR THE LEAST SQUARES METHOD %%

% Specify data to process (from ind_act_error_avg in Main3)
working_data = error_avg_E002;

% Find indices where working_data is not zero and below a threshold
valid_indices = working_data ~= 0 & ~isnan(working_data);% & working_data <10;

% Reshape sigma_array_full and delta_array_full to match working_data
[sigma_array_full, delta_array_full] = meshgrid(sigma_m_array,delta_array);
sigma_array_reshaped = sigma_array_full(valid_indices);
delta_array_reshaped = delta_array_full(valid_indices);
working_data_reshaped = working_data(valid_indices);

% Create a meshgrid from the reshaped data
[sigma_mesh, delta_mesh] = meshgrid(unique(sigma_array_reshaped), unique(delta_array_reshaped));

% Interpolate the error data onto the meshgrid
error_interp = griddata(sigma_array_reshaped, delta_array_reshaped, working_data_reshaped, sigma_mesh, delta_mesh);

% Plot the surface
f = figure;
surf(sigma_mesh, delta_mesh, error_interp, 'FaceColor', 'interp')
%plot3(sigma_array_full(valid_indices),delta_array_full(valid_indices),working_data(valid_indices),'.','Color','blue');

set(gca, 'XScale', 'log','FontSize',textsize)
set(gca, 'YScale', 'log','FontSize',textsize)
set(gca, 'ZScale', 'log','FontSize',textsize)
xlabel('Sigma','FontSize',textsize)
ylabel('Delta','FontSize',textsize)
zlabel('Average Error','FontSize',textsize)
%zlim([0.08 0.09])
title('E0 = 0.2')
colorbar

%% Exponential Error Logic %%
%{
logic_array_exp = zeros(length(delta_array),length(sigma_m_array));
for i = 1:length(delta_array)
    for j = 1:length(sigma_m_array)
        if error_avg_E002(i,j) < 0.305 && error_avg_E002(i,j) > 0 ...
                && error_avg_E003(i,j) < 0.3 && error_avg_E003(i,j) > 0 ...
                && error_avg_E004(i,j) < 0.16 && error_avg_E004(i,j) > 0 ...
                && error_avg_E005(i,j) < 0.1 && error_avg_E005(i,j) > 0 ...
                && error_avg_E01(i,j)  < 0.1 && error_avg_E01(i,j) > 0 ...
                && error_avg_E015(i,j) < 0.1 && error_avg_E015(i,j) > 0 ...
                && error_avg_E02(i,j)  < 0.05 && error_avg_E02(i,j) > 0
            logic_array_exp(i,j) = 1;
        end
    end
end
clear index_exp
[index_exp(:,1),index_exp(:,2)] = find(logic_array_exp==1);
%}

%% Bump-on-tail Error Logic %%
%
logic_array_BOT = zeros(length(delta_array),length(sigma_m_array));
for i = 1:length(delta_array)
    for j = 1:length(sigma_m_array)
        if error_avg_BOT1(i,j) < 0.2 && error_avg_BOT1(i,j) > 0 ...
                && error_avg_BOT2(i,j) < 1 && error_avg_BOT2(i,j) > 0
%                && error_avg_Lin(i,j) < 0.3 && error_avg_Lin(i,j) > 0 ...
            logic_array_BOT(i,j) = 1;
        end
    end
end
clear index_BOT
[index_BOT(:,1),index_BOT(:,2)] = find(logic_array_BOT==1);
%

%% Power Law Error Logic %%
%{
logic_array_POW = zeros(length(delta_array),length(sigma_m_array));
for i = 1:length(delta_array)
    for j = 1:length(sigma_m_array)
        if error_avg_POW_2(i,j) < 0.4 && error_avg_POW_2(i,j) > 0 ...
                && error_avg_POW_3(i,j) < 1 && error_avg_POW_3(i,j) > 0 ...
                && error_avg_POW_4(i,j) < 3 && error_avg_POW_4(i,j) > 0 %...
                && error_avg_POW_5(i,j) < 6   && error_avg_POW_5(i,j) > 0 %...
%                && error_avg_POW_6(i,j) < 150   && error_avg_POW_6(i,j) > 0
            logic_array_POW(i,j) = 1;
        end
    end
end
clear index_POW
[index_POW(:,1),index_POW(:,2)] = find(logic_array_POW==1);
%}

%% Combined Error Logic %%
%{
logic_array_comb = zeros(length(delta_array),length(sigma_m_array));
for i = 1:length(delta_array)
    for j = 1:length(sigma_m_array)
        if logic_array_exp(i,j) == 1 && logic_array_BOT(i,j) == 1
            logic_array_comb(i,j) = 1;
        end
    end
end
clear index_comb
[index_comb(:,1),index_comb(:,2)] = find(logic_array_comb==1);
%}
%{
sig_edges = [0.5:length(sigma_m_array)+0.5];
del_edges = [0.5:length(delta_array)+0.5];
[sig_counts,~] = histcounts(index_POW(:,2),sig_edges);
[del_counts,~] = histcounts(index_POW(:,1),del_edges);


%% Plotting histograms from low-error indices %%
f=figure;
plot(sigma_m_array,sig_counts)
set(gca, 'XScale', 'log','FontSize',textsize)
xlabel('Sigma','FontSize',textsize)

f=figure;
plot(delta_array,del_counts)
set(gca, 'XScale', 'log','FontSize',textsize)
xlabel('Delta','FontSize',textsize)
%}