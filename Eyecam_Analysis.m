%% Paths to stuff
addpath(genpath(cd));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\npy-matlab'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\widefield'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\AP_scripts_cortexlab'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\Lilrig'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\PupilDetection_DLC'));


%% Load DLC output

% beginning of training
% eyecam_DLC_fn = 'C:\Users\Andrada\Desktop\Andy_Pupil DLC\AP108_2021-12-08_2\eyeDLC_resnet50_PupilDetectorApr28shuffle1_1030000.csv';

% end of training
eyecam_DLC_fn = 'C:\Users\Andrada\Desktop\Andy_Pupil DLC\AP108_2021-12-16_2\eyeDLC_resnet50_PupilDetectorApr28shuffle1_1030000.csv';

eyecam_DLC_output = readmatrix(eyecam_DLC_fn);

pupil_top = eyecam_DLC_output(:,2:4); % [x y likelihood]
pupil_bottom = eyecam_DLC_output(:,5:7);
pupil_right = eyecam_DLC_output(:,8:10);
pupil_left = eyecam_DLC_output(:,11:13);
lid_top = eyecam_DLC_output(:,14:16);
lid_bottom = eyecam_DLC_output(:,17:19);
eyecam_DLC_data = [pupil_top, pupil_bottom, pupil_left, pupil_right, lid_top, lid_bottom];

%% Parameters - from pupil DLC page
eyecam_DLC_params.minCertainty = 0.6;
eyecam_DLC_params.lidMinStd = 7;
eyecam_DLC_params.minDistLidCenter = 0.5; % in number of pupil heights
eyecam_DLC_params.surroundingBlinks = 5; % in frames
eyecam_DLC_params.interpMethod = 'linear';
eyecam_DLC_params.maxDistPupilLid = 10; % in pixels
eyecam_DLC_params.smoothSpan = 5; % in frames

%% Some processing steps from pupil DLC page?

% relation between horizontal pupil position, and width and height of pupil
pupil_width = pupil_right(:,1) - pupil_left(:,1);
pupil_height = pupil_bottom(:,2) - pupil_top(:,2);
pupil_centerX = pupil_left(:,1) + 0.5 .* pupil_width;
pupil_valid = all(eyecam_DLC_data(:, 3:3:12) > eyecam_DLC_params.minCertainty, 2);

[F, F0] = dlc.estimateHeightFromWidthPos(pupil_width(pupil_valid), pupil_height(pupil_valid), ...
    pupil_centerX(pupil_valid), eyecam_DLC_params);
    
[pupil_centerY_adj, pupil_height_adj] = dlc.adjustCenterHeight(eyecam_DLC_data, F, eyecam_DLC_params);
    
[blinks, bl_starts, bl_stops] = dlc.detectBlinks(eyecam_DLC_data, eyecam_DLC_params, ...
    pupil_centerY_adj, pupil_height_adj, true);
    
pupil_centerX = medfilt1(pupil_centerX, eyecam_DLC_params.smoothSpan);
pupil_centerY_adj = medfilt1(pupil_centerY_adj, eyecam_DLC_params.smoothSpan);
pupil_center = [pupil_centerX, pupil_centerY_adj];
pupil_center(blinks,:) = NaN;
pupil_height_adj = medfilt1(pupil_height_adj, eyecam_DLC_params.smoothSpan);
diameter = pupil_height_adj;
diameter(blinks) = NaN;

% To do: save data

%% Load exp 

animal = 'AP108';
% day = '2021-12-08';
day = '2021-12-16';
experiment = 2;
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
visual_roi_figure = figure; hold on;
set(gca,'ColorOrder',copper(3));
plot(right_visual_roi.trace');
xlim([1,length(timevec)]);
xline(stim_frame);
legend(num2str(possible_stimuli));

frontal_roi_figure = figure; hold on;
set(gca,'ColorOrder',copper(3));
plot(right_frontal_roi.trace');
xlim([1,length(timevec)]);
xline(stim_frame);
legend(num2str(possible_stimuli));

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

% I think this gives an error because diameter has NaNs in it and it
% probably gets stuck when it tries to interpolate 

% pupil_diameter_stim = interp1(eyecam_t',diameter,time_stimulus);

% use pupil_height_adj instead - still the same error why? 
% - eyecam_t has nans in it too - not sure how to fix (probably to do with experiment start and end time)

% simply find indices of first and last number and use in eyecam_t and
% pupil_height_adj
start_eyecam = find(~isnan(eyecam_t),1,'first');
end_eyecam = find(~isnan(eyecam_t),1,'last');
pupil_diameter_stim = interp1((eyecam_t(start_eyecam:end_eyecam))',pupil_height_adj(start_eyecam:end_eyecam),time_stimulus);

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


figure;hold on;
set(gca,'ColorOrder',copper(3));
for stim_idx =1:length(possible_stimuli)
    these_trials = find(trialStimulusValue==possible_stimuli(stim_idx));
    plot(these_trials, pupil_diameter_stim_avg(these_trials), 'o');
    hold on;
    plot(these_trials, avg_right_visual_roi(these_trials), '*');
    hold on; 
end

figure;hold on;
set(gca,'ColorOrder',copper(3));
for stim_idx =1:length(possible_stimuli)
    plot(stim_idx, pupil_diameter_stim_avg(trialStimulusValue==possible_stimuli(stim_idx)), 'o');
    hold on;
    plot(stim_idx, avg_right_visual_roi(trialStimulusValue==possible_stimuli(stim_idx)), '*');
    hold on; 
end

figure;hold on;
set(gca,'ColorOrder',copper(3));
for stim_idx =1:length(possible_stimuli)
    plot(avg_right_visual_roi(trialStimulusValue==possible_stimuli(stim_idx)), pupil_diameter_stim_avg(trialStimulusValue==possible_stimuli(stim_idx)), 'o');
    hold on; 
end

% plot for each
figure; hold on; 
subplot_idx = 0; 
for stim_idx =1:length(possible_stimuli)
    subplot_idx = subplot_idx + 1;
    subplot (2, 3, subplot_idx)
    plot(avg_right_visual_roi(trialStimulusValue==possible_stimuli(stim_idx)), pupil_diameter_stim_avg(trialStimulusValue==possible_stimuli(stim_idx)), 'o');
    xlabel('Average Visual ROI')
    ylabel('Average pupil diameter')
    title(['Right Visual ROI and stimulus ' num2str(possible_stimuli(stim_idx))])
end

for stim_idx =1:length(possible_stimuli)
    subplot_idx = subplot_idx + 1;
    subplot (2, 3, subplot_idx)
    plot(avg_right_frontal_roi(trialStimulusValue==possible_stimuli(stim_idx)), pupil_diameter_stim_avg(trialStimulusValue==possible_stimuli(stim_idx)), 'o');
    xlabel('Average Visual ROI')
    ylabel('Average pupil diameter')
    title(['Right Frontal ROI and stimulus ' num2str(possible_stimuli(stim_idx))])
end

figure; hold on; 
subplot_idx = 0; 
for stim_idx =1:length(possible_stimuli)
    subplot_idx = subplot_idx + 1;
    subplot (2, 3, subplot_idx)
    plot(avg_left_visual_roi(trialStimulusValue==possible_stimuli(stim_idx)), pupil_diameter_stim_avg(trialStimulusValue==possible_stimuli(stim_idx)), 'o');
    xlabel('Average Visual ROI')
    ylabel('Average pupil diameter')
    title(['Left Visual ROI and stimulus ' num2str(possible_stimuli(stim_idx))])
end

for stim_idx =1:length(possible_stimuli)
    subplot_idx = subplot_idx + 1;
    subplot (2, 3, subplot_idx)
    plot(avg_left_frontal_roi(trialStimulusValue==possible_stimuli(stim_idx)), pupil_diameter_stim_avg(trialStimulusValue==possible_stimuli(stim_idx)), 'o');
    xlabel('Average Visual ROI')
    ylabel('Average pupil diameter')
    title(['Left Frontal ROI and stimulus ' num2str(possible_stimuli(stim_idx))])
end


