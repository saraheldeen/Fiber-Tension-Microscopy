%% Calculate MSD %%
% project: cell-free fluctuation 
% Author: Sarah Eldeen
% Date: Jan/2025
% Goal: Calculate MSDs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Path to files and read the file
clc;clear;close all
%path ='F:\Sarah\FluctuationCellFree\Galvano_Scanner\20240418_noCL_floppy\across_100umabove';
addpath('/Volumes/Sarah_UCI/BEAMS/Fluctuation_Project/20240907_CL_noCL/CL/fibers_find_in_first_5slices_inRef_25slices/Fibers');
addpath('/Volumes/Sarah_UCI/BEAMS/Fluctuation_Project/20240907_CL_noCL/NO_CL_sample2/midway_measurments/F');
addpath('/Volumes/Sarah_UCI/BEAMS/Fluctuation_Project/20240825_fluctuations_NoCL_Sample1');
addpath('/Volumes/Sarah_UCI/BEAMS/Fluctuation_Project/20240822_2rec_far');
addpath('/Volumes/GoogleDrive/My Drive/BEAMS/Matlab_CodesFunctions')
% addpath('D:\Fluctuation_Project\20231101_fluctuation_crosslink_along_across');
after1 = load('T1T2Data.mat'); % in microns 
after2 = load('centroid_T3Z1.mat');
before1 = load('centroid_d50_before_20240907.mat');
before2 = load('centroid_before_20240825.mat');
before3 = load('centroid_far_20240822.mat');
% after2 = load('centroid_T3Z1_updated.mat');
%fter3 = load('kappa.mat');
% Measure the average length of fibers in Fiji. Image already calibrated to
% be in microns (Command + m) = measure 
 % length of fiber in microns
 
 %%
% after20 = load('all_files_centroid_um_gof_only_after20.mat');
% before = load('all_files_centroid_um_gof_only_before.mat');

%%
%after20 = load('files1-8.mat');
%%
figure;

% Initialize plot handles for the legend
h1 = []; % For "After"
h2 = []; % For "Before"
h3 = []; % For "fit"

% Plot "After" data
for i = 1:1%length(after1.all_data)
    centroid_um = after1.all_data{1,i}-after1.all_data{1,i}(1,1);
    lineSpeed =  after1.line_speed{1,1};
    [deltat_test, msdx_test] = calculateMSD(centroid_um, lineSpeed, 1); % in um squared
    msdx_nm2 = msdx_test * 10^6; % in nm squared
    c = polyfit(deltat_test, msdx_nm2, 1);
    y_est = polyval(c, deltat_test);
    kb = 1.381 * 10^-23; % J/K
    T = 296.15; % 23 degree c
    k_spring = 2 * kb * T / c(1); % in Joules/nm^2
    uN_per_um_value = convert_Jnm2_to_uNperum(k_spring);
    
    % Plot and store handles
    h1 = loglog(deltat_test, msdx_nm2, 'k', 'LineWidth', 2);
    hold on;
    h3 = loglog(deltat_test, y_est, 'r--', 'LineWidth', 2);
    spring_constants_cl20(1,i) = uN_per_um_value;
end
%%
% Plot "Before" data
for i = 1:length(before.all_data)
    centroid_um = before.all_data{1,i}(2,:);
    lineSpeed =  (before.all_data{1,i}(1,2)-before.all_data{1,i}(1,1));
    [deltat_test, msdx_test] = calculateMSD(centroid_um, lineSpeed, 1); % in um squared
    msdx_nm2 = msdx_test * 10^6; % in nm squared
    c = polyfit(deltat_test, msdx_nm2, 1);
    y_est = polyval(c, deltat_test);
    kb = 1.381 * 10^-23; % J/K
    T = 296.15; % 23 degree c
    k_spring = 2 * kb * T / c(1); % in Joules/nm^2
    uN_per_um_value = convert_Jnm2_to_uNperum(k_spring);
    
    % Plot and store handles
    h2 = loglog(deltat_test, msdx_nm2, 'b', 'LineWidth', 2);
    loglog(deltat_test, y_est, 'r--', 'LineWidth', 2);
    spring_constants_before(1,i) = uN_per_um_value;
end

ylim([1 10E+3]);
grid on;
xlim([0 1]);

ylabel('$MSD_x$ [nm$^2$]', 'Interpreter', 'latex');
xlabel('$\tau$ [s]', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
ax = gca;
ax.FontSize = 20;

% Add legend
legend([h2(1), h1(1), h3(1)], {'Before', 'After', 'Fit'}, 'Interpreter', 'latex', 'Location', 'best');
%%


%%
%  for file in 20231101_fluctuations_crosslink_along_across
dis_in_um = [39.9 16.4 0 36.1 22.3 43.9 0 56.7];

figure
plot(dis_in_um(1:3), spring_constants_cl20(1:3),'or')
hold on
plot(dis_in_um(6: end), spring_constants_cl20(6:end),'ok')
hold on
plot(dis_in_um(4), spring_constants_cl20(4),'ok')


%%
%%
for i = 1:length(addcl20.all_data)
centroid_um = addcl20.all_data{1,i}(2,:);
lineSpeed =  (addcl20.all_data{1,i}(1,2)-addcl20.all_data{1,i}(1,1));
[deltat_test, msdx_test] = calculateMSD(centroid_um,lineSpeed,1); % in um squared
msdx_nm2 = msdx_test*10^6; % in nm squared
c = polyfit(deltat_test,msdx_nm2,0);
y_est = polyval(c,msdx_nm2);
kb = 1.381*10^-23; %J/K
T = 296.15; %23 degree c
k_spring = 2*kb*T/c; % in Jouls/nm^2
uN_per_um_value = convert_Jnm2_to_uNperum(k_spring);
loglog(deltat_test,msdx_nm2,'k',deltat_test,y_est,'r--','LineWidth',2)
hold all
spring_constants_addcl20(1,i) = uN_per_um_value; 
end
ylim([1 10E+3])
%%
