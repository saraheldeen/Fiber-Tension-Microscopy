%% Calculate minimum frequency required to observe fluctuations %%
% project: cell-free fluctuation 
% Author: Sarah Eldeen
% Date: May/22/2024
% Goal:unknown
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Path to files and read the file
% select a single file directory at a time.
clc;clear;close all
%path ='F:\Sarah\FluctuationCellFree\Galvano_Scanner\20240418_noCL_floppy\across_100umabove';
path = '/Volumes/Sarah_UCI/Fluctuation_Project/20240907/CL/fibers_find_in_first_5slices_inRef_25slices/Fibers';
%path = '/Volumes/Sarah_UCI/Fluctuation_Project/20240907/NO_CL_sample2/midway_measurments/F';
% path = 'D:\Fluctuation_Project\20231113_fluctuation_cl_cellfree_along_across_radial\20231113_CL\across';
%path = 'D:\Fluctuation_Project\20231101_fluctuation_crosslink_along_across';
%path = 'D:\Fluctuation_Project\20231127_fluctuation_CL_Dish4\far';
cd(path);
fileList = dir('T*.oir');
numimg = length(fileList);
pix2um = 1;%/3.2181;
%%
for ii = 1:numimg
    %clearvars -except ii fileList numimg path l1 l2 val data_all line_speed int_img referenceIndex maxLength panData
    fileName = fileList(ii).name;
    fullPath = fullfile(path, fileName);
    file=bfopen(fullPath);
    data = [];
    for i=1:length(file{1, 1})
        data =[data; file{1, 1}{2, 1}];
    end
    A = data;
    data=double(data);
    dataparams=file{1,2};
    lineSpeed=double(string(dataparams.get('speedInformation lineSpeed #2')));
    panX = double(string(dataparams.get('acquisitionValue xPan #1')));
    panY = double(string(dataparams.get('acquisitionValue yPan #1')));
    
    [cycles,position] = size(data);
    %data_m = data-mean(data);
    data_all{ii} = data;
    line_speed{ii} = lineSpeed;
    panData(ii,1) = panX;
    panData(ii,2) = panY;
end

%%

for i = 1: length(file{1,1})
    
    
    
    
end
%% Meta Data Section !!!!!!!!!!!!!!
% stuff=split(string(dataparams),'.');
% %clearvars -except data_all line_speed
% 
% % Get all metadata keys and display them
% keys = dataparams.keySet().iterator();
% while keys.hasNext()
%     key = keys.next();
%     disp(key)
% end
% 
% % Get all metadata keys and display them with their values
% keys = dataparams.keySet().iterator();
% while keys.hasNext()
%     key = keys.next();  % Get the key
%     value = dataparams.get(key);  % Get the corresponding value for the key
%     
%     % Display the key and its value
%     fprintf('Key: %s, Value: %s\n', key, string(value));
% end
% 
% double(string(dataparams.get('acquisitionValue xPan #1')))
% double(string(dataparams.get('acquisitionValue yPan #1')))% 
% 
% %axis startPosition #2
% %
% %axis basePosition #4
% X1 = double(string(dataparams.get('axis startPosition #1')))
% Y1 = double(string(dataparams.get('axis endPosition #1')))
% 
% X2 = double(string(dataparams.get('axis startPosition #2')))
% Y2 = double(string(dataparams.get('axis endPosition #2')))
% 
% X3 = double(string(dataparams.get('axis startPosition #3')))
% Y3 = double(string(dataparams.get('axis endPosition #3')))
% 
% X4 = double(string(dataparams.get('axis startPosition #4')))
% Y4 = double(string(dataparams.get('axis endPosition #4')))

%%
%%
figure
for i = 1%:length(data_all)
    imagesc(data_all{1,i}(1:500:end,:),[0 700])
    lags = 1:1:length(data_all{1,i});
    space = 1:1:size(data_all{1,i},2);
end
ylabel('Time [s]', 'Interpreter','latex');
xlabel('Space [$\mu$m]', 'Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
ax = gca;ax.FontSize = 20;%axis equal
max_y_labels = 20;
lags_in_seconds =lags*line_speed{i}/1000;
step_size = ceil(length(lags_in_seconds) / max_y_labels);
label_indices = 1:step_size:length(lags_in_seconds);
yticks(label_indices);
yticklabels(arrayfun(@(i) sprintf('%.1f', lags_in_seconds(label_indices(i))), 1:length(label_indices), 'UniformOutput', false));

space_um = space*pix2um;
step_size = ceil(length(space) / max_y_labels);
label_indices = 1:step_size:length(space);
xticks(label_indices);
xticklabels(arrayfun(@(i) sprintf('%.1f', space_um(label_indices(i))), 1:length(label_indices), 'UniformOutput', false));

%%
% all_data{1,1} = add;
%%  string(dataparams).split(',')
%% This section is for centroid tracking by fitting a gaussian
lim = 7; % is in pixels 8 for before CL
%pix2um = 1 / 32.1812; % pix/um
% Preallocate cell arrays
num_data = length(data_all);
all_goodness_of_fit = cell(1, num_data);
all_data = cell(1, num_data);
% figure
downsample = 50;
%% Updates Gaussian fit with limits of 0.8
for ii = 1:num_data
    lineSpeed = line_speed{1};
    data = data_all{ii}(1:downsample:end,:);
    cycles = size(data, 1);
    
    % Preallocate arrays for centroid and goodness_of_fit
    centroid = zeros(1, cycles);
    goodness_of_fit = zeros(1, cycles);
    
    for i = 1:cycles
        input = data(i, :);
        [maxValue, maxIndex] = max(input);
        ind1 = max(maxIndex - lim, 1);
        ind2 = min(maxIndex + lim, length(input));
        input_segment = input(ind1:ind2);
        x = 1:numel(input_segment);
        
        % Fit Gaussian with bounds to avoid fitting issues
        try
            fitOptions = fitoptions('gauss1', 'Lower', [0, 0, 0], 'Upper', [Inf, numel(input_segment), Inf]);
            [gaussFit, gof] = fit(x', input_segment', 'gauss1', fitOptions);
            coeffs = coeffvalues(gaussFit);
            mu = coeffs(2); % Mean of the Gaussian
            centroid(i) = mu;
            goodness_of_fit(i) = gof.rsquare;
            
            % Check the goodness of fit
            if goodness_of_fit(i) < 0.8
                centroid(i) = NaN; % Replace centroid with NaN if goodness of fit is less than 0.8
            end
            
            % Plotting
%             plot(x, input_segment, '.k', 'MarkerSize', 20); % Plot data points in black with specified marker size
%             hold on;
%             plot(x, gaussFit(x), '-r', 'LineWidth', 3); % Plot fitted curve in red with specified line width
%             legend({'Raw Data', 'Gaussian fit'}, 'Interpreter', 'latex', 'Location', 'NorthEast');
            
        catch ME
            % Handle fit errors (e.g., Inf in model function)
            centroid(i) = NaN; % Assign NaN if fit fails
            goodness_of_fit(i) = NaN;
        end
    end
    
    % Convert centroid to microns
    centroid_um = centroid * pix2um;
    
    % Store results
    all_goodness_of_fit{ii} = goodness_of_fit; % Store goodness of fit for each cycle
    all_data{ii} = [centroid_um]; % Store centroid_um in all_data
    
    disp(ii);
end

% % Labeling and axes formatting
% ylabel('Intensity [a.u.]', 'Interpreter', 'latex');
% xlabel('Space [$\mu$m]', 'Interpreter', 'latex');
% set(gca, 'TickLabelInterpreter', 'latex');
% ax = gca;
% ax.FontSize = 20;
% 
% % X-axis ticks
% max_y_labels = 8;
% space = 1:1:length(input_segment);
% space_um = space * pix2um;
% step_size = ceil(length(space) / max_y_labels);
% label_indices = 1:step_size:length(space);
% xticks(label_indices);
% xticklabels(arrayfun(@(i) sprintf('%.1f', space_um(label_indices(i))), 1:length(label_indices), 'UniformOutput', false));
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Below is extra 
%%
% figure
for ii = 1:num_data
    lineSpeed = line_speed{ii};
    data = data_all{ii}(1:downsample:end,:);
    cycles = size(data, 1);
    % Preallocate arrays for centroid and goodness_of_fit
    centroid = zeros(1, cycles);
    goodness_of_fit = zeros(1, cycles);
    
    for i = 1:cycles
        input = data(i, :);
        [maxValue, maxIndex] = max(input);
        ind1 = max(maxIndex - lim, 1);
        ind2 = min(maxIndex + lim, length(input));
        input_segment = input(ind1:ind2);
        x = 1:numel(input_segment);
        
        % Fit Gaussian with bounds to avoid fitting issues
        try
            fitOptions = fitoptions('gauss1', 'Lower', [0, 0, 0], 'Upper', [Inf, numel(input_segment), Inf]);
            [gaussFit, gof] = fit(x', input_segment', 'gauss1', fitOptions);
            coeffs = coeffvalues(gaussFit);
            mu = coeffs(2); % Mean of the Gaussian
            centroid(i) = mu;
            goodness_of_fit(i) = gof.rsquare;
%             plot(x, input_segment, '.k', 'MarkerSize', 20); % Plot data points in black with specified marker size
%             hold on;
%             plot(x, gaussFit(x), '-r', 'LineWidth', 3); % Plot fitted curve in red with specified line width
%             legend({'Raw Data', 'Gaussian fit'}, 'Interpreter', 'latex', 'Location', 'NorthEast');
        catch ME
            % Handle fit errors (e.g., Inf in model function)
            %disp(['Fit error at cycle ' num2str(i) ': ' ME.message]);
            centroid(i) = NaN; % Assign NaN if fit fails
            goodness_of_fit(i) = NaN;
        end
    end
    
 
        centroid_um = centroid ;%* pix2um;
        all_goodness_of_fit{ii} = goodness_of_fit; % Store goodness of fit for each cycle
        all_data{ii} = [centroid_um]; % Store centroid_um in all_data

    disp(ii);
end

% ylabel('Intensity [a.u.]', 'Interpreter','latex');
% xlabel('Space [$\mu$m]', 'Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% ax = gca;ax.FontSize = 20;
% max_y_labels = 8;
% space = 1:1:length(input_segment);
% space_um = space*pix2um;
% step_size = ceil(length(space) / max_y_labels);
% label_indices = 1:step_size:length(space);
% xticks(label_indices);
% xticklabels(arrayfun(@(i) sprintf('%.1f', space_um(label_indices(i))), 1:length(label_indices), 'UniformOutput', false));

%%
% All data contains centroid in microns!! 
%%

figure 
time = (1:1:length(centroid))*downsample*lineSpeed / 1000; % in seconds
centroid_um = centroid * pix2um;
plot(time, centroid_um)
ylabel('centroid position [um]', 'Interpreter','latex');
xlabel('time [s]', 'Interpreter','latex');
ax = gca;ax.FontSize = 20;
set(gca,'TickLabelInterpreter','latex');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculates MSD of centroid 
[deltat, msdx] = calculateMSD(centroid_um,1,1);

figure
loglog(deltat,msdx,'o')
xlabel('$t$ [s]', 'Interpreter','latex')
ylabel('$MSDx$ $[um^2]$', 'Interpreter','latex')
ax = gca;
ax.FontSize = 20;
grid on
set(gca,'TickLabelInterpreter','latex')
%%
% ylabel('Intensity Value', 'Interpreter','latex');
% xlabel('Pixel', 'Interpreter','latex');
% ax = gca;ax.FontSize = 20;
% set(gca,'TickLabelInterpreter','latex');
% Define the path and extract the last three folders
folders = strsplit(path, '/');
filename = sprintf('%s_%s_%s_%s.mat', folders{end-2}, folders{end-1}, folders{end}, 'data');
% Specify the full file path
full_path = '/Users/saraheldeen/Desktop/';
% Save the variables to the specified path with the generated filename
save(fullfile(full_path, filename), 'all_data');
%%