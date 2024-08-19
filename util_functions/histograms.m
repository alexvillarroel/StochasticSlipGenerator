% specify the folders
tic;
% folders = {'C:\Users\axlph\OneDrive - Universidad de Concepción\magister\Proyecto de Tesis\Slips\Simulation_8.7',...
%     'C:\Users\axlph\OneDrive - Universidad de Concepción\magister\Proyecto de Tesis\Slips\Simulation_8.8',...
%     'C:\Users\axlph\OneDrive - Universidad de Concepción\magister\Proyecto de Tesis\Slips\Simulation_8.9',...
%     'C:\Users\axlph\OneDrive - Universidad de Concepción\magister\Proyecto de Tesis\Slips\Simulation_9.0',...
%     'C:\Users\axlph\OneDrive - Universidad de Concepción\magister\Proyecto de Tesis\Slips\Simulation_9.1',...
%     'C:\Users\axlph\OneDrive - Universidad de Concepción\magister\Proyecto de Tesis\Slips\Simulation_9.2'};
folders = {'C:\Users\axlph\OneDrive - Universidad de Concepción\magister\Proyecto de Tesis\Slips\Simulation_9.3'};
% create a figure
figure;

% iterate over each folder
for i = 1:length(folders)
    % get all .mat files in the folder
    mat_files = dir(fullfile(folders{i}, '*.mat'));

    % initialize an array to store all max slips
    max_slips = [];

    % iterate over each file
    for j = 1:length(mat_files)
        % load the .mat file
        mat_data = load(fullfile(folders{i}, mat_files(j).name));
    
        % assuming 'slip' is a variable in the .mat file
        if isfield(mat_data, 'slip')
            % get the maximum of 'slip'
            max_slip = max(mat_data.slip(:));
            max_slips = [max_slips; max_slip];
        end
    end

    % calculate mean and standard deviation
    mean_slip = median(max_slips);
    std_slip = std(max_slips);

    % create a subplot for this folder
    subplot(2, 3, i);
    histogram(max_slips);
    hold on;

    % plot mean and standard deviations
    mean_line = xline(mean_slip, 'r');
    std1_line = xline(mean_slip + std_slip, 'g');
    xline(mean_slip - std_slip, 'g');
    std2_line = xline(mean_slip + 2*std_slip, 'b');
    xline(mean_slip - 2*std_slip, 'b');
    std3_line = xline(mean_slip + 3*std_slip, 'k');
    xline(mean_slip - 3*std_slip, 'k');

    % set the plot box aspect ratio

    % add a legend
    legend([mean_line, std1_line, std2_line, std3_line], {'Median', 'Median ± 1 STD', 'Median ± 2 STD', 'Median ± 3 STD'}, 'Location', 'best');
    ylim([0,1200])
    xlabel('Max Slip[m]')
    ylabel('Number of distributions')
    % set the plot box aspect ratio
    pbaspect([1 1 1]);

    hold off;
end
elapsed_time=toc;
fprintf('elapsed time is %.2f seconds.\n',elapsed_time)