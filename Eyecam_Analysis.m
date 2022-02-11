%% Paths to stuff
addpath(genpath(cd));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\npy-matlab'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\widefield'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\Lilrig'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\PupilDetection_DLC'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\AP_scripts_cortexlab'));

%% Load exp 

animal = 'AP108';
% day = '2021-12-08';
day = '2021-12-16';
experiment = 2;
load_parts.imaging = true;
load_parts.cam = true;
verbose = true;
AP_load_experiment;

%% Get activity for 2 seconds after stim onset in ROIs 

% create matrix of times for movie
timestep = 0.01;
start_time = 0;
end_time = 2;
timevec = [start_time:timestep:end_time];

stim_frame = (-start_time)*(1/timestep)+1;

time_stimulus = stimOn_times+timevec;

% find activity for above time
all_stim_act = interp1(frame_t,fVdf',time_stimulus);
all_stim_act = permute(all_stim_act, [3,2,1]);
all_stim_act = all_stim_act - all_stim_act(:,stim_frame,:);

% get fluoresence
all_stim_interval_fluorescence = AP_svdFrameReconstruct(Udf,all_stim_act);

% right visual ROI
AP_image_scroll(all_stim_interval_fluorescence);
axis image;

right_visual_roi = roi;

% right frontal ROI
AP_image_scroll(all_stim_interval_fluorescence,timevec);
axis image;

right_frontal_roi = roi;

% left visual ROI
AP_image_scroll(all_stim_interval_fluorescence);
axis image;

left_visual_roi = roi;

% right frontal ROI
AP_image_scroll(all_stim_interval_fluorescence,timevec);
axis image;

left_frontal_roi = roi;

% plot to check - legend doesn't work! 
% visual_roi_figure = figure; hold on;
% set(gca,'ColorOrder',copper(3));
% plot(right_visual_roi.trace');
% xlim([1,length(timevec)]);
% xline(stim_frame);
% legend(num2str(possible_stimuli));
% 
% frontal_roi_figure = figure; hold on;
% set(gca,'ColorOrder',copper(3));
% plot(right_frontal_roi.trace');
% xlim([1,length(timevec)]);
% xline(stim_frame);
% legend(num2str(possible_stimuli));

% average activity 

avg_right_visual_roi = nanmean(right_visual_roi.trace, 2);
avg_right_frontal_roi = nanmean(right_frontal_roi.trace,2);
avg_left_visual_roi = nanmean(left_visual_roi.trace, 2);
avg_left_frontal_roi = nanmean(left_frontal_roi.trace,2);

% plot check
% figure; hold on;
% plot(avg_visual_roi,'bo')
% hold on;
% plot(avg_frontal_roi,'ro')

%% Get average pupil diameter from eyecam data

start_eyecam = find(~isnan(eyecam_t),1,'first');
end_eyecam = find(~isnan(eyecam_t),1,'last');
pupil_diameter_stim = interp1((eyecam_t(start_eyecam:end_eyecam))',pupil_diameter(start_eyecam:end_eyecam),time_stimulus);

% do average 
pupil_diameter_stim_avg = nanmean(pupil_diameter_stim,2);

%% plot of averages per trial 

% find possible stimuli
possible_contrasts = unique(signals_events.stimContrastValues);
possible_sides = unique(signals_events.stimAzimuthValues)/90;
possible_stimuli = possible_sides.*possible_contrasts;
possible_stimuli = unique(possible_stimuli);

% find value of stimulus per trial
trialStimulusValue = signals_events.stimAzimuthValues/90 .* signals_events.stimContrastValues;


% figure;hold on;
% set(gca,'ColorOrder',copper(3));
% for stim_idx =1:length(possible_stimuli)
%     these_trials = find(trialStimulusValue==possible_stimuli(stim_idx));
%     plot(these_trials, pupil_diameter_stim_avg(these_trials), 'o');
%     hold on;
%     plot(these_trials, avg_right_visual_roi(these_trials), '*');
%     hold on; 
% end
% 
% figure;hold on;
% set(gca,'ColorOrder',copper(3));
% for stim_idx =1:length(possible_stimuli)
%     plot(stim_idx, pupil_diameter_stim_avg(trialStimulusValue==possible_stimuli(stim_idx)), 'o');
%     hold on;
%     plot(stim_idx, avg_right_visual_roi(trialStimulusValue==possible_stimuli(stim_idx)), '*');
%     hold on; 
% end
% 
% figure;hold on;
% set(gca,'ColorOrder',copper(3));
% for stim_idx =1:length(possible_stimuli)
%     plot(avg_right_visual_roi(trialStimulusValue==possible_stimuli(stim_idx)), pupil_diameter_stim_avg(trialStimulusValue==possible_stimuli(stim_idx)), 'o');
%     hold on; 
% end
% 
% % plot for each
% figure; hold on; 
% subplot_idx = 0; 
% for stim_idx =1:length(possible_stimuli)
%     subplot_idx = subplot_idx + 1;
%     subplot (2, 3, subplot_idx)
%     plot(avg_right_visual_roi(trialStimulusValue==possible_stimuli(stim_idx)), pupil_diameter_stim_avg(trialStimulusValue==possible_stimuli(stim_idx)), 'o');
%     xlabel('Average Visual ROI')
%     ylabel('Average pupil diameter')
%     title(['Right Visual ROI and stimulus ' num2str(possible_stimuli(stim_idx))])
% end
% 
% for stim_idx =1:length(possible_stimuli)
%     subplot_idx = subplot_idx + 1;
%     subplot (2, 3, subplot_idx)
%     plot(avg_right_frontal_roi(trialStimulusValue==possible_stimuli(stim_idx)), pupil_diameter_stim_avg(trialStimulusValue==possible_stimuli(stim_idx)), 'o');
%     xlabel('Average Visual ROI')
%     ylabel('Average pupil diameter')
%     title(['Right Frontal ROI and stimulus ' num2str(possible_stimuli(stim_idx))])
% end
% 
% figure; hold on; 
% subplot_idx = 0; 
% for stim_idx =1:length(possible_stimuli)
%     subplot_idx = subplot_idx + 1;
%     subplot (2, 3, subplot_idx)
%     plot(avg_left_visual_roi(trialStimulusValue==possible_stimuli(stim_idx)), pupil_diameter_stim_avg(trialStimulusValue==possible_stimuli(stim_idx)), 'o');
%     xlabel('Average Visual ROI')
%     ylabel('Average pupil diameter')
%     title(['Left Visual ROI and stimulus ' num2str(possible_stimuli(stim_idx))])
% end
% 
% for stim_idx =1:length(possible_stimuli)
%     subplot_idx = subplot_idx + 1;
%     subplot (2, 3, subplot_idx)
%     plot(avg_left_frontal_roi(trialStimulusValue==possible_stimuli(stim_idx)), pupil_diameter_stim_avg(trialStimulusValue==possible_stimuli(stim_idx)), 'o');
%     xlabel('Average Visual ROI')
%     ylabel('Average pupil diameter')
%     title(['Left Frontal ROI and stimulus ' num2str(possible_stimuli(stim_idx))])
% end

%% 3 plots
% one fluorescence, one pupil,  - mean
% and line in between stuff

figure;hold on;
set(gca,'ColorOrder',copper(1));
pupil_diameter_stim_mean = nan(length(possible_stimuli),1);
for stim_idx =1:length(possible_stimuli)
    plot(stim_idx, pupil_diameter_stim_avg(trialStimulusValue==possible_stimuli(stim_idx)), 'o');
    hold on;
    pupil_diameter_stim_mean(stim_idx) = mean(pupil_diameter_stim_avg(trialStimulusValue==possible_stimuli(stim_idx)));
end
% add mean to plot
plot(pupil_diameter_stim_mean, 's', 'MarkerSize', 15, 'MarkerFaceColor', 'k')
plot(pupil_diameter_stim_mean)
title('Pupil diameter')

% p-value is 0.3
ranksum(pupil_diameter_stim_avg(trialStimulusValue==-1),pupil_diameter_stim_avg(trialStimulusValue==1))


% Check
% figure
% boxplot(pupil_diameter_stim_avg', trialStimulusValue, 'PlotStyle','compact'); 
% hold on;
% plot(pupil_diameter_stim_mean);

% Right Visual
figure;hold on;
set(gca,'ColorOrder',copper(3));
avg_right_visual_roi_stim_mean = nan(length(possible_stimuli),1);
for stim_idx =1:length(possible_stimuli)
    plot(stim_idx, avg_right_visual_roi(trialStimulusValue==possible_stimuli(stim_idx)), 'o');
    avg_right_visual_roi_stim_mean(stim_idx) = mean(avg_right_visual_roi(trialStimulusValue==possible_stimuli(stim_idx)));
end
% add mean to plot
plot(avg_right_visual_roi_stim_mean, 's', 'MarkerSize', 15, 'MarkerFaceColor', 'k')
plot(avg_right_visual_roi_stim_mean)
title('Right Visual')

% if window of time was smaller this would make more sense
ranksum(avg_right_visual_roi(trialStimulusValue==-1),avg_right_visual_roi(trialStimulusValue==1))

% Left Visual
figure;hold on;
set(gca,'ColorOrder',copper(3));
avg_left_visual_roi_stim_mean = nan(length(possible_stimuli),1);
for stim_idx =1:length(possible_stimuli)
    plot(stim_idx, avg_left_visual_roi(trialStimulusValue==possible_stimuli(stim_idx)), 'o');
    avg_left_visual_roi_stim_mean(stim_idx) = mean(avg_left_visual_roi(trialStimulusValue==possible_stimuli(stim_idx)));
end
% add mean to plot
plot(avg_left_visual_roi_stim_mean, 's', 'MarkerSize', 15, 'MarkerFaceColor', 'k')
plot(avg_left_visual_roi_stim_mean)
title('Left Visual')

% again smaller window!!
ranksum(avg_left_visual_roi(trialStimulusValue==-1),avg_left_visual_roi(trialStimulusValue==1))


% Right Frontal
figure;hold on;
set(gca,'ColorOrder',copper(3));
avg_right_frontal_roi_stim_mean = nan(length(possible_stimuli),1);
for stim_idx =1:length(possible_stimuli)
    plot(stim_idx, avg_right_frontal_roi(trialStimulusValue==possible_stimuli(stim_idx)), 'o');
    avg_right_frontal_roi_stim_mean(stim_idx) = mean(avg_right_frontal_roi(trialStimulusValue==possible_stimuli(stim_idx)));
end
% add mean to plot
plot(avg_right_frontal_roi_stim_mean, 's', 'MarkerSize', 15, 'MarkerFaceColor', 'k')
plot(avg_right_frontal_roi_stim_mean)
title('Right Frontal')

% p is 0.51
ranksum(avg_right_frontal_roi(trialStimulusValue==-1),avg_right_frontal_roi(trialStimulusValue==1))


% Left Frontal
figure;hold on;
set(gca,'ColorOrder',copper(3));
avg_left_frontal_roi_stim_mean = nan(length(possible_stimuli),1);
for stim_idx =1:length(possible_stimuli)
    plot(stim_idx, avg_left_frontal_roi(trialStimulusValue==possible_stimuli(stim_idx)), 'o');
    avg_left_frontal_roi_stim_mean(stim_idx) = mean(avg_left_frontal_roi(trialStimulusValue==possible_stimuli(stim_idx)));
end
% add mean to plot
plot(avg_left_frontal_roi_stim_mean, 's', 'MarkerSize', 15, 'MarkerFaceColor', 'k')
plot(avg_left_frontal_roi_stim_mean)
title('Left Frontal')

% p is 0.29
ranksum(avg_left_frontal_roi(trialStimulusValue==-1),avg_left_frontal_roi(trialStimulusValue==1))


% HOW DO I GET IT TO DO EACH STIM IN A DIFFERENT COLOUR?

% combined (z-score before) - use errorbar instead of dots: standard dev or
% standard error of mean (std/sqrt(#points))

figure;hold on;
for stim_idx =1:length(possible_stimuli)
    avg_right_visual_roi_stim = avg_right_visual_roi(trialStimulusValue==possible_stimuli(stim_idx));
    [avg_right_visual_roi_stim_norm, avg_right_visual_roi_stim_mu, avg_right_visual_roi_stim_sigma] = zscore(avg_right_visual_roi_stim);
    avg_right_visual_roi_stim_std_error_mean = avg_right_visual_roi_stim_sigma/sqrt(length(avg_right_visual_roi_stim));
    plot(stim_idx, avg_right_visual_roi_stim_norm, 'o', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
    hold on;
%     plot(stim_idx, avg_right_visual_roi_stim_mu, 's', 'MarkerSize', 15, 'MarkerFaceColor', 'k');
%     hold on;
end

for stim_idx =1:length(possible_stimuli)
    pupil_diameter_stim = pupil_diameter_stim_avg(trialStimulusValue==possible_stimuli(stim_idx));
    [pupil_diameter_stim_norm, pupil_diameter_stim_mu, pupil_diameter_stim_sigma] = zscore(pupil_diameter_stim);
    pupil_diameter_stim_std_error_mean = pupil_diameter_stim_sigma/sqrt(length(pupil_diameter_stim));
    plot(stim_idx, pupil_diameter_stim_norm, 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
%     errorbar(stim_idx, pupil_diameter_stim_norm, pupil_diameter_stim_std_error_mean, ...
%         '-s', 'MarkerSize',10, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
    hold on;
%     plot(stim_idx, pupil_diameter_stim_mu, 's', 'MarkerSize', 15, 'MarkerFaceColor', 'k');
%     hold on;
end

% errorbar(x,y,err,'-s','MarkerSize',10,...
%     'MarkerEdgeColor','red','MarkerFaceColor','red')



%% Other plot
% like this but label with stuff from correlation and do all areas
% for correlation
% corrcoef - get R-value and p-value

% columns are stimuli and rows are areas - sgtitle('Subplot Grid Title') ?
figure;hold on;
subplot_idx = 0;

for stim_idx =1:length(possible_stimuli)
    subplot_idx = subplot_idx + 1;
    subplot (4, 3, subplot_idx)
    plot(avg_right_visual_roi(trialStimulusValue==possible_stimuli(stim_idx)), pupil_diameter_stim_avg(trialStimulusValue==possible_stimuli(stim_idx)), 'o');
    [R, p] = corrcoef(avg_right_visual_roi(trialStimulusValue==possible_stimuli(stim_idx)), pupil_diameter_stim_avg(trialStimulusValue==possible_stimuli(stim_idx)));
    hold on; 
    title(['R = ', num2str(R(1,2)), ' and p = ', num2str(p(1,2))])
end

for stim_idx =1:length(possible_stimuli)
    subplot_idx = subplot_idx + 1;
    subplot (4, 3, subplot_idx)
    plot(avg_left_visual_roi(trialStimulusValue==possible_stimuli(stim_idx)), pupil_diameter_stim_avg(trialStimulusValue==possible_stimuli(stim_idx)), 'o');
    [R, p] = corrcoef(avg_left_visual_roi(trialStimulusValue==possible_stimuli(stim_idx)), pupil_diameter_stim_avg(trialStimulusValue==possible_stimuli(stim_idx)));
    hold on; 
    title(['R = ', num2str(R(1,2)), ' and p = ', num2str(p(1,2))])
end

for stim_idx =1:length(possible_stimuli)
    subplot_idx = subplot_idx + 1;
    subplot (4, 3, subplot_idx)
    plot(avg_right_frontal_roi(trialStimulusValue==possible_stimuli(stim_idx)), pupil_diameter_stim_avg(trialStimulusValue==possible_stimuli(stim_idx)), 'o');
    [R, p] = corrcoef(avg_right_frontal_roi(trialStimulusValue==possible_stimuli(stim_idx)), pupil_diameter_stim_avg(trialStimulusValue==possible_stimuli(stim_idx)));
    hold on; 
    title(['R = ', num2str(R(1,2)), ' and p = ', num2str(p(1,2))])
end

for stim_idx =1:length(possible_stimuli)
    subplot_idx = subplot_idx + 1;
    subplot (4, 3, subplot_idx)
    plot(avg_left_frontal_roi(trialStimulusValue==possible_stimuli(stim_idx)), pupil_diameter_stim_avg(trialStimulusValue==possible_stimuli(stim_idx)), 'o');
    [R, p] = corrcoef(avg_left_frontal_roi(trialStimulusValue==possible_stimuli(stim_idx)), pupil_diameter_stim_avg(trialStimulusValue==possible_stimuli(stim_idx)));
    hold on; 
    title(['R = ', num2str(R(1,2)), ' and p = ', num2str(p(1,2))])
end

%% New plot 
% across days - mean pupil for each mouse per learning day

