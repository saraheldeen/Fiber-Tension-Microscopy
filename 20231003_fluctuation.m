%% this code is to analyze point scans and check how downsampling affects the parameter to be observed

clc;
clear all;
cd('/Volumes/Sarah_UCI/Fluctuation_Project/pointscans_CL_20231127')
addpath('/Volumes/GoogleDrive/My Drive/BEAMS/Matlab_CodesFunctions')  


T = table2array(readtable('tile9_FOV_z5_X120Y120_LS10fibers_alongacross_5PS_0021.csv')); % QPD Signal 
time = T(:,2);
signal = T(:,3)-mean(T(:,3));

% Downsampling section 
y1 = signal(1:500000,1);
t1 = time(1:500000,1);%1:100000

y2 = downsample(y1,2);
t2 = downsample(t1,2);

y3 = downsample(y1,3);
t3 = downsample(t1,3);

y4 = downsample(y1,4);
t4 = downsample(t1,4);

y5 = downsample(y1,5);
t5 = downsample(t1,5);

y10 = downsample(y1,10);
t10 = downsample(t1,10);

y100 = downsample(y1,100);
t100 = downsample(t1,100);

y500 = downsample(y1,500);
t500 = downsample(t1,500);

y1000 = downsample(y1,1000);
t1000 = downsample(t1,1000);

y2000 = downsample(y1,2000);
t2000 = downsample(t1,2000);
%%
% plotting
figure
plot(t1/1000,y1,'ok')
hold all
% plot(t2/1000,y2,'.b')
% plot(t3/1000,y3,'.r')
% plot(t4/1000,y4,'.g')
% plot(t5/1000,y5,'.y')
 plot(t10/1000,y10,'.c')
% plot(t100/1000,y100,'.')
% plot(t500/1000,y500,'.')
% plot(t1000/1000,y1000,'.')
plot(t2000/1000,y2000,'*r')
legend('raw data', 'DS 10', 'DS 2000','Interpreter','latex');
xlabel('Lags [s]', 'Interpreter','latex');
ylabel('Intensity', 'Interpreter','latex');
%%
% Do not run it will take forever
% [R1, TimeLag1] = CalculateXCorr(y1, y1, t1);
% [R2, TimeLag2] = CalculateXCorr(y2, y2, t2);
% [R3, TimeLag3] = CalculateXCorr(y3, y3, t3);
% [R4, TimeLag4] = CalculateXCorr(y4, y4, t4);
% [R5, TimeLag5] = CalculateXCorr(y5, y5, t5);
% [R10, TimeLag10] = CalculateXCorr(y10, y10, t10);
% [R100, TimeLag100] = CalculateXCorr(y100, y100, t100);
% [R500, TimeLag500] = CalculateXCorr(y500, y500, t500);
% [R1000, TimeLag1000] = CalculateXCorr(y1000, y1000, t1000);
% [R2500, TimeLag2500] = CalculateXCorr(y2500, y2500, t2500);
% using xcorr 
[tR1, TimeLag1] = xcorr(y1, y1,'Normalized');
non_negative_lag_indices = find(TimeLag1 >= 0);
R1= tR1(non_negative_lag_indices);
lags1 = TimeLag1(non_negative_lag_indices);

[tR2, TimeLag2] = xcorr(y2, y2,'Normalized');
non_negative_lag_indices = find(TimeLag2 >= 0);
R2= tR2(non_negative_lag_indices);
lags2 = TimeLag2(non_negative_lag_indices);

[tR3, TimeLag3] = xcorr(y3, y3,'Normalized');
non_negative_lag_indices = find(TimeLag3 >= 0);
R3= tR3(non_negative_lag_indices);
lags3 = TimeLag3(non_negative_lag_indices);

[tR4, TimeLag4] = xcorr(y4, y4,'Normalized');
non_negative_lag_indices = find(TimeLag4 >= 0);
R4= tR4(non_negative_lag_indices);
lags4 = TimeLag4(non_negative_lag_indices);

[tR5, TimeLag5] = xcorr(y5, y5,'Normalized');
non_negative_lag_indices = find(TimeLag5 >= 0);
R5= tR5(non_negative_lag_indices);
lags5 = TimeLag5(non_negative_lag_indices);

[tR10, TimeLag10] = xcorr(y10, y10,'Normalized');
non_negative_lag_indices = find(TimeLag10 >= 0);
R10= tR10(non_negative_lag_indices);
lags10 = TimeLag10(non_negative_lag_indices);

[tR100, TimeLag100] = xcorr(y100, y100,'Normalized');
non_negative_lag_indices = find(TimeLag100 >= 0);
R100= tR100(non_negative_lag_indices);
lags100 = TimeLag100(non_negative_lag_indices);

[tR500, TimeLag500] = xcorr(y500, y500,'Normalized');
non_negative_lag_indices = find(TimeLag500 >= 0);
R500= tR500(non_negative_lag_indices);
lags500 = TimeLag500(non_negative_lag_indices);

[tR1000, TimeLag1000] = xcorr(y1000, y1000,'Normalized');
non_negative_lag_indices = find(TimeLag1000 >= 0);
R1000= tR1000(non_negative_lag_indices);
lags1000 = TimeLag1000(non_negative_lag_indices);

[tR2000, TimeLag2000] = xcorr(y2000, y2000,'Normalized');
non_negative_lag_indices = find(TimeLag2000 >= 0);
R2000= tR2000(non_negative_lag_indices);
lags2000 = TimeLag2000(non_negative_lag_indices);
%
% figure
% cf = 2/1000;
% plot(lags1000*cf,R1000,'.-b')
% hold on 
% plot(lags1000*cf,R1000s,'*r')
%% 
figure
plot(t1/1000,log(R1+1),'.-')
hold all
% plot(t2/1000,log(R2+1),'.-')
% plot(t4/1000,log(R4+1),'.-')
% plot(t5/1000,log(R5+1),'.-')
% plot(t10/1000,log(R10+1),'.-')
plot(t100/1000,log(R500+1),'.-')
% plot(t1000/1000,log(R1000+1),'.-')
plot(t2000/1000,log(R2000+1),'.-')
%
legend('fps 500,000','fps 1000','fps 250','Interpreter','latex');
 
%legend('fps 500,000','fps 250,000','fps 125,000','fps 100,000','fps 50,000','fps 5,000','fps 500','fps 250','Interpreter','latex');
 xlabel('Lags [s]', 'Interpreter','latex');
 ylabel('log(corr+1)', 'Interpreter','latex');

%%
% figure
% plot(t1/1000,(R1),'.-')
% hold all
% plot(t2/1000,(R2),'.-')
% plot(t4/1000,(R4),'.-')
% plot(t5/1000,(R5),'.-')
% plot(t10/1000,(R10),'.-')
% plot(t100/1000,(R100),'.-')
% plot(t500/1000,(R500),'.-')
% plot(t1000/1000,(R1000),'.-')
% plot(t2000/1000,(R2000),'.-')
% 
% legend('fps 500,000','fps 250,000','fps 125,000','fps 100,000','fps 50,000','fps 5,000','fps 1,000','fps 500','fps 250','Interpreter','latex');
% xlabel('Lags [s]', 'Interpreter','latex');
% ylabel('Auto corr', 'Interpreter','latex');

%
l1 = -0.3; l2 = 0.3;max_y_labels = 4;
font = 18;
figure;
rows = 9;
subplot(rows,1,1)
imagesc((R1),[l1 l2]);
title('fps 500,000', 'Interpreter','latex');
lags_in_seconds = t1/1000;

ylabel('Lags [s]', 'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex')
colormap jet;
c = colorbar;
c.TickLabelInterpreter = 'latex';
ax = gca;ax.FontSize = font; 
set(gca,'XTick',[])
step_size = ceil(length(lags_in_seconds) / max_y_labels);label_indices = 1:step_size:length(lags_in_seconds);
yticks(label_indices);
yticklabels(arrayfun(@(i) sprintf('%.2f', lags_in_seconds(label_indices(i))), 1:length(label_indices), 'UniformOutput', false));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(rows,1,2)
imagesc((R2),[l1 l2]);
title('fps 250,000', 'Interpreter','latex');
lags_in_seconds = t2/1000;

ylabel('Lags [s]', 'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex')
colormap jet;
c = colorbar;
c.TickLabelInterpreter = 'latex';
set(gca,'XTick',[])
ax = gca;ax.FontSize = font; 
step_size = ceil(length(lags_in_seconds) / max_y_labels);label_indices = 1:step_size:length(lags_in_seconds);
yticks(label_indices);
yticklabels(arrayfun(@(i) sprintf('%.2f', lags_in_seconds(label_indices(i))), 1:length(label_indices), 'UniformOutput', false));
ax.XAxis.TickLength = [0 0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(rows,1,3)
imagesc((R4),[l1 l2]);
title('fps 125,000', 'Interpreter','latex');
lags_in_seconds = t4/1000;

ylabel('Lags [s]', 'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex')
colormap jet;
c = colorbar;
c.TickLabelInterpreter = 'latex';
ax = gca;ax.FontSize = font; 
step_size = ceil(length(lags_in_seconds) / max_y_labels);label_indices = 1:step_size:length(lags_in_seconds);
yticks(label_indices);
yticklabels(arrayfun(@(i) sprintf('%.2f', lags_in_seconds(label_indices(i))), 1:length(label_indices), 'UniformOutput', false));
set(gca,'XTick',[])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(rows,1,4)
imagesc((R5),[l1 l2]);
title('fps 100,000', 'Interpreter','latex');
lags_in_seconds = t5/1000;

ylabel('Lags [s]', 'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex')
colormap jet;
c = colorbar;
c.TickLabelInterpreter = 'latex';
ax = gca;ax.FontSize = font; 
step_size = ceil(length(lags_in_seconds) / max_y_labels);label_indices = 1:step_size:length(lags_in_seconds);
yticks(label_indices);
yticklabels(arrayfun(@(i) sprintf('%.2f', lags_in_seconds(label_indices(i))), 1:length(label_indices), 'UniformOutput', false));
set(gca,'XTick',[])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(rows,1,5)
imagesc((R10),[l1 l2]);
title('fps 50,000', 'Interpreter','latex');
lags_in_seconds = t10/1000;

ylabel('Lags [s]', 'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex')
colormap jet;
c = colorbar;
c.TickLabelInterpreter = 'latex';
ax = gca;ax.FontSize = font; 
step_size = ceil(length(lags_in_seconds) / max_y_labels);label_indices = 1:step_size:length(lags_in_seconds);
yticks(label_indices);
yticklabels(arrayfun(@(i) sprintf('%.2f', lags_in_seconds(label_indices(i))), 1:length(label_indices), 'UniformOutput', false));
set(gca,'XTick',[])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(rows,1,6)
imagesc((R100),[l1 l2]);
title('fps 5,000', 'Interpreter','latex');
lags_in_seconds = t100/1000;

ylabel('Lags [s]', 'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex')
colormap jet;
c = colorbar;
c.TickLabelInterpreter = 'latex';
ax = gca;ax.FontSize = font; 
step_size = ceil(length(lags_in_seconds) / max_y_labels);label_indices = 1:step_size:length(lags_in_seconds);
yticks(label_indices);
yticklabels(arrayfun(@(i) sprintf('%.2f', lags_in_seconds(label_indices(i))), 1:length(label_indices), 'UniformOutput', false));
set(gca,'XTick',[])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(rows,1,7)
imagesc((R500),[l1 l2]);
title('fps 1,000', 'Interpreter','latex');
lags_in_seconds = t500/1000;

ylabel('Lags [s]', 'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex')
colormap jet;
c = colorbar;
c.TickLabelInterpreter = 'latex';
ax = gca;ax.FontSize = font; 
step_size = ceil(length(lags_in_seconds) / max_y_labels);label_indices = 1:step_size:length(lags_in_seconds);
yticks(label_indices);
yticklabels(arrayfun(@(i) sprintf('%.2f', lags_in_seconds(label_indices(i))), 1:length(label_indices), 'UniformOutput', false));
set(gca,'XTick',[])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(rows,1,8)
imagesc((R1000),[l1 l2]);
title('fps 500', 'Interpreter','latex');
lags_in_seconds = t1000/ 1000;

ylabel('Lags [s]', 'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex')
colormap jet;
c = colorbar;
c.TickLabelInterpreter = 'latex';
ax = gca;ax.FontSize = font; 
step_size = ceil(length(lags_in_seconds) / max_y_labels);label_indices = 1:step_size:length(lags_in_seconds);
yticks(label_indices);
yticklabels(arrayfun(@(i) sprintf('%.2f', lags_in_seconds(label_indices(i))), 1:length(label_indices), 'UniformOutput', false));
set(gca,'XTick',[])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(rows,1,9)
imagesc((R2000),[l1 l2]);
title('fps 250', 'Interpreter','latex');
lags_in_seconds = t2000/ 1000;

ylabel('Lags [s]', 'Interpreter','latex');
set(gca,'TickLabelInterpreter','latex')
colormap jet;
c = colorbar;
c.TickLabelInterpreter = 'latex';
ax = gca;ax.FontSize = font; 
step_size = ceil(length(lags_in_seconds) / max_y_labels);label_indices = 1:step_size:length(lags_in_seconds);
yticks(label_indices);
yticklabels(arrayfun(@(i) sprintf('%.2f', lags_in_seconds(label_indices(i))), 1:length(label_indices), 'UniformOutput', false));
set(gca,'XTick',[])


%%
% imagesc(CC_pos,[l1 l2]);axis image
% ylabel('Auto Correlation at Lags (seconds)', 'Interpreter','latex');xlabel('Space', 'Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% colormap jet;axis equal
% c = colorbar;c.TickLabelInterpreter = 'latex';
% ax = gca;ax.FontSize = 20;daspect([5 1000 1])%axis equal
% max_y_labels = 10;lags_in_seconds = lags * lineSpeed / 1000;
% step_size = ceil(length(lags_in_seconds) / max_y_labels);label_indices = 1:step_size:length(lags_in_seconds);
% yticks(label_indices);
% yticklabels(arrayfun(@(i) sprintf('%.1f', lags_in_seconds(label_indices(i))), 1:length(label_indices), 'UniformOutput', false));

%%
% l1 = -0.3; l2 = 0.3;max_y_labels = 3;
% figure;
% rows = 8
% subplot(rows,1,1)
% imagesc(log(R1+1),[l1 l2]);
% title('fps 500,000', 'Interpreter','latex');
% lags_in_seconds = t1/1000;
% 
% ylabel('Lags [s]', 'Interpreter','latex');
% set(gca,'TickLabelInterpreter','latex')
% colormap jet;
% c = colorbar;
% c.TickLabelInterpreter = 'latex';
% ax = gca;ax.FontSize = 20; 
% set(gca,'XTick',[])
% step_size = ceil(length(lags_in_seconds) / max_y_labels);label_indices = 1:step_size:length(lags_in_seconds);
% yticks(label_indices);
% yticklabels(arrayfun(@(i) sprintf('%.2f', lags_in_seconds(label_indices(i))), 1:length(label_indices), 'UniformOutput', false));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(rows,1,2)
% imagesc(log(R2+1),[l1 l2]);
% title('fps 250,000', 'Interpreter','latex');
% lags_in_seconds = t2/1000;
% 
% ylabel('Lags [s]', 'Interpreter','latex');
% set(gca,'TickLabelInterpreter','latex')
% colormap jet;
% c = colorbar;
% c.TickLabelInterpreter = 'latex';
% set(gca,'XTick',[])
% ax = gca;ax.FontSize = 20; 
% step_size = ceil(length(lags_in_seconds) / max_y_labels);label_indices = 1:step_size:length(lags_in_seconds);
% yticks(label_indices);
% yticklabels(arrayfun(@(i) sprintf('%.2f', lags_in_seconds(label_indices(i))), 1:length(label_indices), 'UniformOutput', false));
% ax.XAxis.TickLength = [0 0];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(rows,1,3)
% imagesc(log(R4+1),[l1 l2]);
% title('fps 125,000', 'Interpreter','latex');
% lags_in_seconds = t4/1000;
% 
% ylabel('Lags [s]', 'Interpreter','latex');
% set(gca,'TickLabelInterpreter','latex')
% colormap jet;
% c = colorbar;
% c.TickLabelInterpreter = 'latex';
% ax = gca;ax.FontSize = 20; 
% step_size = ceil(length(lags_in_seconds) / max_y_labels);label_indices = 1:step_size:length(lags_in_seconds);
% yticks(label_indices);
% yticklabels(arrayfun(@(i) sprintf('%.2f', lags_in_seconds(label_indices(i))), 1:length(label_indices), 'UniformOutput', false));
% set(gca,'XTick',[])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(rows,1,4)
% imagesc(log(R5+1),[l1 l2]);
% title('fps 100,000', 'Interpreter','latex');
% lags_in_seconds = t5/1000;
% 
% ylabel('Lags [s]', 'Interpreter','latex');
% set(gca,'TickLabelInterpreter','latex')
% colormap jet;
% c = colorbar;
% c.TickLabelInterpreter = 'latex';
% ax = gca;ax.FontSize = 20; 
% step_size = ceil(length(lags_in_seconds) / max_y_labels);label_indices = 1:step_size:length(lags_in_seconds);
% yticks(label_indices);
% yticklabels(arrayfun(@(i) sprintf('%.2f', lags_in_seconds(label_indices(i))), 1:length(label_indices), 'UniformOutput', false));
% set(gca,'XTick',[])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(rows,1,5)
% imagesc(log(R10+1),[l1 l2]);
% title('fps 50,000', 'Interpreter','latex');
% lags_in_seconds = t10/1000;
% 
% ylabel('Lags [s]', 'Interpreter','latex');
% set(gca,'TickLabelInterpreter','latex')
% colormap jet;
% c = colorbar;
% c.TickLabelInterpreter = 'latex';
% ax = gca;ax.FontSize = 20; 
% step_size = ceil(length(lags_in_seconds) / max_y_labels);label_indices = 1:step_size:length(lags_in_seconds);
% yticks(label_indices);
% yticklabels(arrayfun(@(i) sprintf('%.2f', lags_in_seconds(label_indices(i))), 1:length(label_indices), 'UniformOutput', false));
% set(gca,'XTick',[])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(rows,1,6)
% imagesc(log(R100+1),[l1 l2]);
% title('fps 5,000', 'Interpreter','latex');
% lags_in_seconds = t100/1000;
% 
% ylabel('Lags [s]', 'Interpreter','latex');
% set(gca,'TickLabelInterpreter','latex')
% colormap jet;
% c = colorbar;
% c.TickLabelInterpreter = 'latex';
% ax = gca;ax.FontSize = 20; 
% step_size = ceil(length(lags_in_seconds) / max_y_labels);label_indices = 1:step_size:length(lags_in_seconds);
% yticks(label_indices);
% yticklabels(arrayfun(@(i) sprintf('%.2f', lags_in_seconds(label_indices(i))), 1:length(label_indices), 'UniformOutput', false));
% set(gca,'XTick',[])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(rows,1,7)
% imagesc(log(R1000+1),[l1 l2]);
% title('fps 500', 'Interpreter','latex');
% lags_in_seconds = t1000/ 1000;
% 
% ylabel('Lags [s]', 'Interpreter','latex');
% set(gca,'TickLabelInterpreter','latex')
% colormap jet;
% c = colorbar;
% c.TickLabelInterpreter = 'latex';
% ax = gca;ax.FontSize = 20; 
% step_size = ceil(length(lags_in_seconds) / max_y_labels);label_indices = 1:step_size:length(lags_in_seconds);
% yticks(label_indices);
% yticklabels(arrayfun(@(i) sprintf('%.2f', lags_in_seconds(label_indices(i))), 1:length(label_indices), 'UniformOutput', false));
% set(gca,'XTick',[])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(rows,1,8)
% imagesc(log(R2000+1),[l1 l2]);
% title('fps 250', 'Interpreter','latex');
% lags_in_seconds = t2000/ 1000;
% 
% ylabel('Lags [s]', 'Interpreter','latex');
% set(gca,'TickLabelInterpreter','latex')
% colormap jet;
% c = colorbar;
% c.TickLabelInterpreter = 'latex';
% ax = gca;ax.FontSize = 20; 
% step_size = ceil(length(lags_in_seconds) / max_y_labels);label_indices = 1:step_size:length(lags_in_seconds);
% yticks(label_indices);
% yticklabels(arrayfun(@(i) sprintf('%.2f', lags_in_seconds(label_indices(i))), 1:length(label_indices), 'UniformOutput', false));
% set(gca,'XTick',[])
% 

%%
% imagesc(CC_pos,[l1 l2]);axis image
% ylabel('Auto Correlation at Lags (seconds)', 'Interpreter','latex');xlabel('Space', 'Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% colormap jet;axis equal
% c = colorbar;c.TickLabelInterpreter = 'latex';
% ax = gca;ax.FontSize = 20;daspect([5 1000 1])%axis equal
% max_y_labels = 10;lags_in_seconds = lags * lineSpeed / 1000;
% step_size = ceil(length(lags_in_seconds) / max_y_labels);label_indices = 1:step_size:length(lags_in_seconds);
% yticks(label_indices);
% yticklabels(arrayfun(@(i) sprintf('%.1f', lags_in_seconds(label_indices(i))), 1:length(label_indices), 'UniformOutput', false));

%%
% %%
% figure;
% rows = 2
% subplot(rows,1,1)
% imagesc(log(R1+1));
% colormap('hot'); % You can change the colormap to your preference
% colorbar; % Display a colorbar
% axis off; % Turn off axis labels and ticks
% title('full');
%axis([1 length(R1) 1 2]); % You can adjust the aspect ratio by changing the axis limits

%%
% % Create a heatmap
% figure;
% rows = 9
% subplot(rows,1,1)
% imagesc(log(y1));
% colormap('hot'); % You can change the colormap to your preference
% colorbar; % Display a colorbar
% axis off; % Turn off axis labels and ticks
% title('full');
% axis([1 length(y1) 1 2]); % You can adjust the aspect ratio by changing the axis limits
% 
% 
% subplot(rows,1,2)
% imagesc(y2);
% colormap('hot'); % You can change the colormap to your preference
% colorbar; % Display a colorbar
% axis off; % Turn off axis labels and ticks
% title('downsampled by 2');
% % Adjust the aspect ratio to make it look like a straight line
% axis([1 length(y2) 1 2]); % You can adjust the aspect ratio by changing the axis limits
% 
% subplot(rows,1,3)
% imagesc(y3);
% colormap('hot'); % You can change the colormap to your preference
% colorbar; % Display a colorbar
% axis off; % Turn off axis labels and ticks
% title('downsampled by 3');
% % Adjust the aspect ratio to make it look like a straight line
% axis([1 length(y3) 1 2]); % You can adjust the aspect ratio by changing the axis limits
% 
% subplot(rows,1,4)
% imagesc(y4);
% colormap('hot'); % You can change the colormap to your preference
% colorbar; % Display a colorbar
% axis off; % Turn off axis labels and ticks
% title('downsampled by 4');
% % Adjust the aspect ratio to make it look like a straight line
% axis([1 length(y4) 1 2]); % You can adjust the aspect ratio by changing the axis limits
% 
% subplot(rows,1,5)
% imagesc(y5);
% colormap('hot'); % You can change the colormap to your preference
% colorbar; % Display a colorbar
% axis off; % Turn off axis labels and ticks
% title('downsampled by 5');
% % Adjust the aspect ratio to make it look like a straight line
% axis([1 length(y5) 1 2]); % You can adjust the aspect ratio by changing the axis limits
% 
% subplot(rows,1,6)
% imagesc(y10);
% colormap('hot'); % You can change the colormap to your preference
% colorbar; % Display a colorbar
% axis off; % Turn off axis labels and ticks
% title('downsampled by 10');
% % Adjust the aspect ratio to make it look like a straight line
% axis([1 length(y10) 1 2]); % You can adjust the aspect ratio by changing the axis limits
% 
% subplot(rows,1,7)
% imagesc(y100);
% colormap('hot'); % You can change the colormap to your preference
% colorbar; % Display a colorbar
% axis off; % Turn off axis labels and ticks
% title('downsampled by 100');
% % Adjust the aspect ratio to make it look like a straight line
% axis([1 length(y100) 1 2]); % You can adjust the aspect ratio by changing the axis limits
% 
% subplot(rows,1,8)
% imagesc(y500);
% colormap('hot'); % You can change the colormap to your preference
% colorbar; % Display a colorbar
% axis off; % Turn off axis labels and ticks
% title('downsampled by 500');
% % Adjust the aspect ratio to make it look like a straight line
% axis([1 length(y500) 1 2]); % You can adjust the aspect ratio by changing the axis limits
% 
% subplot(rows,1,9)
% imagesc(y1000);
% colormap('hot'); % You can change the colormap to your preference
% colorbar; % Display a colorbar
% axis off; % Turn off axis labels and ticks
% title('downsampled by 1000');
% % Adjust the aspect ratio to make it look like a straight line
% axis([1 length(y1000) 1 2]); % You can adjust the aspect ratio by changing the axis limits
% 
% %%
% figure
% [R1, TimeLag1] = CalculateXCorr(signal, signal, time)
% plot(TimeLag1,R1)
% %%
% figure
% plot(signal)
% %%
% figure
% plot(TimeLag1,R1)
% hold all
% plot(TimeLag2,R2)
% plot(TimeLag3,R3)
% plot(TimeLag4,R4)
% plot(TimeLag5,R5)
% plot(TimeLag10,R10)
% legend
% xlim([-50 50])


