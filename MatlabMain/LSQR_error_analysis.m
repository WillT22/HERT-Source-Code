working_data = ind_act_error_avg;

% Find indices where working_data is not zero
valid_indices = working_data ~= 0 & ~isnan(working_data) & working_data <1;

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

set(gca, 'XScale', 'log','FontSize',textsize)
set(gca, 'YScale', 'log','FontSize',textsize)
set(gca, 'ZScale', 'log','FontSize',textsize)
xlabel('Sigma','FontSize',textsize)
ylabel('Delta','FontSize',textsize)
zlabel('Average Error','FontSize',textsize)
%zlim([0.08 0.09])
title('BOT1')
colorbar

%{
logic_array_exp = zeros(length(delta_array),length(sigma_m_array));
max_value = 1.0;
for i = 1:length(delta_array)
    for j = 1:length(sigma_m_array)
        if error_avg_E002(i,j) < 0.8 && error_avg_E002(i,j) > 0 ...
                && error_avg_E003(i,j) < 0.4 && error_avg_E003(i,j) > 0 ...
                && error_avg_E004(i,j) < 1 && error_avg_E004(i,j) > 0 ...
                && error_avg_E005(i,j) < 0.11 && error_avg_E005(i,j) > 0 ...
                && error_avg_E01(i,j)  < 0.2 && error_avg_E01(i,j) > 0 ...
                && error_avg_E015(i,j) < 0.1 && error_avg_E015(i,j) > 0 ...
                && error_avg_E02(i,j) < 0.13 && error_avg_E02(i,j) > 0 ...
            logic_array_exp(i,j) = 1;
        end
    end
end
clear index_exp
[index_exp(:,1),index_exp(:,2)] = find(logic_array_exp==1);
%}
logic_array_BOT = zeros(length(delta_array),length(sigma_m_array));
for i = 1:length(delta_array)
    for j = 1:length(sigma_m_array)
        if error_avg_BOT1(i,j) < 0.18 && error_avg_BOT1(i,j) > 0 ...
                && error_avg_BOT2(i,j) < 0.3 && error_avg_BOT2(i,j) > 0
%                && error_avg_Lin(i,j) < 0.3 && error_avg_Lin(i,j) > 0 ...
            logic_array_BOT(i,j) = 1;
        end
    end
end
clear index_BOT
[index_BOT(:,1),index_BOT(:,2)] = find(logic_array_BOT==1);
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
[sig_counts,~] = histcounts(index_exp(:,2),sig_edges);
[del_counts,~] = histcounts(index_exp(:,1),del_edges);

f=figure;
plot(sigma_m_array,sig_counts)
set(gca, 'XScale', 'log','FontSize',textsize)

f=figure;
plot(delta_array,del_counts)
set(gca, 'XScale', 'log','FontSize',textsize)
%}