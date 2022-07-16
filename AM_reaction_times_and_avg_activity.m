%% Load and plot reaction times

%% all mice
load('all_mice_behaviour.mat');

animals = {'AP107','AP108', 'AP113', 'AP114', 'AP115'};

% Original task plot
figure;
for animal_id=1:length(animals)
    all_training_days = behaviour(animal_id).day;
    original_task_days_mask = behaviour(animal_id).original_task_days_mask;
    original_task_day_idx = find(original_task_days_mask);
    % get reaction times percentage that falls between 100-200 ms
    reaction_times_percentage = nan(1,length(original_task_day_idx));
    first_day_idx = original_task_day_idx(1);
    for day_idx = original_task_day_idx
        reaction_times = behaviour(animal_id).reaction_times{day_idx};
        reaction_times_percentage(day_idx-first_day_idx+1) = sum(discretize(reaction_times,[0.1 0.2])==1)/length(reaction_times);
    end
    % make plot
    plot(reaction_times_percentage,'o-')
    hold on;
end

%title('Reaction times for original task')
xlabel('Training day')
ylabel('% of fast reaction times')
legend(animals, 'Location','northwest')
AM_prettyfig

% Reversal task plot
figure;
for animal_id=1:length(animals)
    all_training_days = behaviour(animal_id).day;
    reversal_task_days_mask = behaviour(animal_id).reversal_task_days_mask;
    reversal_task_day_idx = find(reversal_task_days_mask);
    first_day_idx = reversal_task_day_idx(1);
    % get reaction times percentage that falls between 100-200 ms
    reaction_times_percentage = nan(1,length(reversal_task_day_idx));
    for day_idx = reversal_task_day_idx
        reaction_times = behaviour(animal_id).reaction_times{day_idx};
        reaction_times_percentage(day_idx-first_day_idx+1) = sum(discretize(reaction_times,[0.1 0.2])==1)/length(reaction_times);
    end
    % make plot
    plot(1:length(reversal_task_day_idx),reaction_times_percentage,'o-')
    hold on;
end

%title('Reaction times for reversal task')
xlabel('Training day')
ylabel('% of fast reaction times')
legend(animals, 'Location','northwest')
AM_prettyfig

%% AP107 histograms
load('all_mice_behaviour.mat');

animals = {'AP107','AP108', 'AP113', 'AP114', 'AP115'};

animal_id = 1;

% original
learned_day_idx = 5;
% get day indexes
all_training_days = behaviour(animal_id).day;
original_task_days_mask = behaviour(animal_id).original_task_days_mask;
original_task_day_idx = find(original_task_days_mask);
first_day_idx = original_task_day_idx(1);
% find day indices for learned days
these_day_idx = original_task_day_idx(learned_day_idx:end);
% append all RT post-learning
original_RT = [];
for day_idx = these_day_idx
    reaction_times = behaviour(animal_id).reaction_times{day_idx};
    original_RT = [original_RT reaction_times];
end

% reversal
learned_day_idx = 7;
% get day indexes
all_training_days = behaviour(animal_id).day;
reversal_task_days_mask = behaviour(animal_id).reversal_task_days_mask;
reversal_task_day_idx = find(reversal_task_days_mask);
first_day_idx = reversal_task_day_idx(1);
% find day indices for learned days
these_day_idx = reversal_task_day_idx(learned_day_idx:end);
% append all RT post-learning
reversal_RT = [];
for day_idx = these_day_idx
    reaction_times = behaviour(animal_id).reaction_times{day_idx};
    reversal_RT = [reversal_RT reaction_times];
end

% plot
figure; 
histogram(original_RT,'BinWidth', 0.01, 'BinLimits', [0 1], 'Normalization', 'probability')
hold on;
histogram(reversal_RT,'BinWidth', 0.01, 'BinLimits', [0 1], 'Normalization', 'probability')
legend({'right stim / left move task', 'left stim / left move task'})
hold on;
xlabel('Reaction time')
ylabel('Probability')
AM_prettyfig
%% PASSIVE FLUORESCENCE ----------------------------------------------------

%% ROIs all mice
% load activity from passive
load('passive_master_activity.mat');
animals = {passive_master_activity.animal};

timestep = passive_master_activity.timestep;
start_time = passive_master_activity.start_time;
end_time = passive_master_activity.end_time;
stim_frame = passive_master_activity.stim_frame;
timevec = passive_master_activity.timevec;
% num_comp = passive_master_activity.num_comp;
num_comp = 200;

% load ROI masks
load('ROIs.mat');

% load master U
master_u_fn = 'D:\Andrada\Master project\widefield_alignment\U_master';
load(master_u_fn);


%% - frontal left ROIs 
% get day indexes for 107
animal_id = 1;
all_training_days = passive_master_activity(animal_id).day;
original_task_days_mask = passive_master_activity(animal_id).original_task_days_mask;
muscimol_days_mask = passive_master_activity(animal_id).muscimol_days_mask;
original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);

frontal_left_avg_ROI = nan(length(original_task_day_idx),length(timevec),length(animals));
for animal_id=1:length(animals)
    all_training_days = passive_master_activity(animal_id).day;
    original_task_days_mask = passive_master_activity(animal_id).original_task_days_mask;
    muscimol_days_mask = passive_master_activity(animal_id).muscimol_days_mask;
    original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);
    
    for day_idx=original_task_day_idx
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
        deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
        
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        right_stim_act = deconvolved_stim_avg_act(:,:,trialStimulusValue==1);
        frontal_left = AP_svd_roi(U_master(:,:,1:num_comp),right_stim_act,[],[],roi.frontal_left);
        frontal_left = permute(frontal_left,[2,3,1]);
        frontal_left_avg_ROI(day_idx-original_task_day_idx(1)+1,:,animal_id) = mean(frontal_left,2);
    end
end

% average across animals
avg_frontal_left_ROI = nanmean(frontal_left_avg_ROI,3);
std_frontal_left = std(frontal_left_avg_ROI, 0, 3);
sem_frontal_left =  std_frontal_left./sqrt(length(animals));

% make a plot of just day 1, learned day 5 and last day 10
days = [1 5 10];
roi_colors = max(0, brewermap(3,'Reds')-0.2);
figure('Name', 'Passive original frontal left ROI');
hold on;
set(gca,'ColorOrder',roi_colors)
AP_errorfill(timevec,avg_frontal_left_ROI(days,:)',sem_frontal_left(days,:)',roi_colors)
ylabel('Fluorescence')
xlabel('Time from stimulus onset')
legend({'', 'day 1', '', 'day 5', '', 'day 10'})
AM_prettyfig

%% - frontal right ROIs
% get day indexes for 115
animal_id = 5;
all_training_days = passive_master_activity(animal_id).day;
reversal_task_days_mask = passive_master_activity(animal_id).reversal_task_days_mask;
reversal_task_day_idx = find(reversal_task_days_mask);

frontal_right_avg_ROI = nan(length(reversal_task_day_idx),length(timevec),length(animals));
for animal_id=1:length(animals)
    
    all_training_days = passive_master_activity(animal_id).day;
    reversal_task_days_mask = passive_master_activity(animal_id).reversal_task_days_mask;
    reversal_task_day_idx = find(reversal_task_days_mask);
    
    for day_idx=reversal_task_day_idx
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
        deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
        
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        right_stim_act = deconvolved_stim_avg_act(:,:,trialStimulusValue==-1);
        frontal_right = AP_svd_roi(U_master(:,:,1:num_comp),right_stim_act,[],[],roi.frontal_right);
        frontal_right = permute(frontal_right,[2,3,1]);
        frontal_right_avg_ROI(day_idx-reversal_task_day_idx(1)+1,:,animal_id) = mean(frontal_right,2);
    end
end

% average across non-learning animals
avg_frontal_right_ROI = nanmean(frontal_right_avg_ROI(:,:,2:end),3);
std_frontal_right = std(frontal_right_avg_ROI(:,:,2:end), 0, 3);
sem_frontal_right =  std_frontal_right./sqrt(length(animals));

% make a plot of just day 1, learned day 5 and last day 10
days = [1 5 10];
roi_colors = max(0, brewermap(3,'Blues')-0.2);
figure('Name', 'Passive original frontal right ROI');
hold on;
set(gca,'ColorOrder',roi_colors)
AP_errorfill(timevec,avg_frontal_right_ROI(days,:)',sem_frontal_right(days,:)',roi_colors)
ylabel('Fluorescence')
xlabel('Time from stimulus onset')
legend({'', 'day 1', '', 'day 5', '', 'day 10'})
ylim([-2.5*10^-3 2.5*10^-3])
AM_prettyfig

% make a plot for 107
days = [1 5 10];
roi_colors = max(0, brewermap(3,'Blues')-0.2);
figure('Name', 'Passive reversal frontal right ROI learner');
hold on;
set(gca,'ColorOrder',roi_colors)
plot(timevec,frontal_right_avg_ROI(days, :, 1)', 'LineWidth', 1.5);
ylabel('Fluorescence')
xlabel('Time from stimulus onset')
ylim([-2.5*10^-3 2.5*10^-3])
legend({'day 1', 'day 5', 'day 10'})
AM_prettyfig

%% - visual left ROIs 
% get day indexes for 107
animal_id = 1;
all_training_days = passive_master_activity(animal_id).day;
original_task_days_mask = passive_master_activity(animal_id).original_task_days_mask;
muscimol_days_mask = passive_master_activity(animal_id).muscimol_days_mask;
original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);

visual_left_avg_ROI = nan(length(original_task_day_idx),length(timevec),length(animals));
for animal_id=1:length(animals)
    all_training_days = passive_master_activity(animal_id).day;
    original_task_days_mask = passive_master_activity(animal_id).original_task_days_mask;
    muscimol_days_mask = passive_master_activity(animal_id).muscimol_days_mask;
    original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);
    
    for day_idx=original_task_day_idx
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
        deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
        
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        right_stim_act = deconvolved_stim_avg_act(:,:,trialStimulusValue==1);
        visual_left = AP_svd_roi(U_master(:,:,1:num_comp),right_stim_act,[],[],roi.visual_left);
        visual_left = permute(visual_left,[2,3,1]);
        visual_left_avg_ROI(day_idx-original_task_day_idx(1)+1,:,animal_id) = mean(visual_left,2);
    end
end

% average across animals
avg_visual_left_ROI = nanmean(visual_left_avg_ROI,3);
std_visual_left = std(visual_left_avg_ROI, 0, 3);
sem_visual_left =  std_visual_left./sqrt(length(animals));

% make a plot of just day 1, learned day 5 and last day 10
days = [1 5 10];
roi_colors = brewermap(3,'Oranges');
figure('Name', 'Passive original visual left ROI');
hold on;
set(gca,'ColorOrder',roi_colors)
AP_errorfill(timevec,avg_visual_left_ROI(days,:)',sem_visual_left(days,:)',roi_colors)
ylabel('Fluorescence')
xlabel('Time from stimulus onset')
legend({'', 'day 1', '', 'day 5', '', 'day 10'})
AM_prettyfig

%% - visual right ROIs
% get day indexes for 115
animal_id = 5;
all_training_days = passive_master_activity(animal_id).day;
reversal_task_days_mask = passive_master_activity(animal_id).reversal_task_days_mask;
reversal_task_day_idx = find(reversal_task_days_mask);

visual_right_avg_ROI = nan(length(reversal_task_day_idx),length(timevec),length(animals));
for animal_id=1:length(animals)
    
    all_training_days = passive_master_activity(animal_id).day;
    reversal_task_days_mask = passive_master_activity(animal_id).reversal_task_days_mask;
    reversal_task_day_idx = find(reversal_task_days_mask);
    
    for day_idx=reversal_task_day_idx
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
        deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
        
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        right_stim_act = deconvolved_stim_avg_act(:,:,trialStimulusValue==-1);
        visual_right = AP_svd_roi(U_master(:,:,1:num_comp),right_stim_act,[],[],roi.visual_right);
        visual_right = permute(visual_right,[2,3,1]);
        visual_right_avg_ROI(day_idx-reversal_task_day_idx(1)+1,:,animal_id) = mean(visual_right,2);
    end
end

% average across non-learning animals
avg_visual_right_ROI = nanmean(visual_right_avg_ROI(:,:,2:end),3);
std_visual_right = std(visual_right_avg_ROI(:,:,2:end), 0, 3);
sem_visual_right =  std_visual_right./sqrt(length(animals));

% make a plot of just day 1, learned day 5 and last day 10
days = [1 5 10];
roi_colors = max(0, brewermap(3,'Purples')-0.2);
figure('Name', 'Passive original visual right ROI');
hold on;
set(gca,'ColorOrder',roi_colors)
AP_errorfill(timevec,avg_visual_right_ROI(days,:)',sem_visual_right(days,:)',roi_colors)
ylabel('Fluorescence')
xlabel('Time from stimulus onset')
ylim([-5*10^-3 5*10^-3])
legend({'', 'day 1', '', 'day 5', '', 'day 10'})
AM_prettyfig

% make plot for 107
days = [1 5 10];
roi_colors = max(0, brewermap(3,'Purples')-0.2);
figure('Name', 'Passive reversal visual right ROI learner');
hold on;
set(gca,'ColorOrder',roi_colors)
plot(timevec,visual_right_avg_ROI(days, :, 1)', 'LineWidth', 1.5);
ylabel('Fluorescence')
xlabel('Time from stimulus onset')
ylim([-5*10^-3 5*10^-3])
legend({'day 1', 'day 5', 'day 10'})
AM_prettyfig


%% - post-training vis and front

% - visual left ROIs 
% get day indexes for 113
animal_id = 3;
all_training_days = passive_master_activity(animal_id).day;
post_training_passive_days_mask = passive_master_activity(animal_id).post_training_passive_days_mask;
post_training_passive_day_idx = find(post_training_passive_days_mask);

visual_left_avg_ROI = nan(length(post_training_passive_day_idx),length(timevec),length(animals));
for animal_id=3:length(animals)
    all_training_days = passive_master_activity(animal_id).day;
    post_training_passive_days_mask = passive_master_activity(animal_id).post_training_passive_days_mask;
    post_training_passive_day_idx = find(post_training_passive_days_mask);

    for day_idx=post_training_passive_day_idx
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
        deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
        
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        right_stim_act = deconvolved_stim_avg_act(:,:,trialStimulusValue==1);
        visual_left = AP_svd_roi(U_master(:,:,1:num_comp),right_stim_act,[],[],roi.visual_left);
        visual_left = permute(visual_left,[2,3,1]);
        visual_left_avg_ROI(day_idx-post_training_passive_day_idx(1)+1,:,animal_id) = mean(visual_left,2);
    end
end

% average across animals
this_visual_left_avg_ROI = visual_left_avg_ROI(:,:,3:end);
avg_visual_left_ROI = nanmean(this_visual_left_avg_ROI,3);
std_visual_left = std(this_visual_left_avg_ROI, 0, 3);
sem_visual_left =  std_visual_left./sqrt(length(animals)-2);


% get last day of training ROIs

last_day_training_idx = post_training_passive_day_idx(1)-1;

last_visual_left_ROI = nan(length(timevec),length(animals));
day_idx = last_day_training_idx;
for animal_id=3:length(animals)
    stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
    deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
    deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
    
    trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
    right_stim_act = deconvolved_stim_avg_act(:,:,trialStimulusValue==1);
    visual_left = AP_svd_roi(U_master(:,:,1:num_comp),right_stim_act,[],[],roi.visual_left);
    visual_left = permute(visual_left,[2,3,1]);
    last_visual_left_ROI(:,animal_id) = mean(visual_left,2);
end

% average across animals
this_visual_left_avg_ROI = last_visual_left_ROI(:,3:end);
avg_last_visual_left_ROI = nanmean(this_visual_left_avg_ROI,2);
std_last_visual_left = std(this_visual_left_avg_ROI, 0, 2);
sem_last_visual_left =  std_last_visual_left./sqrt(length(animals)-2);

% make a plot of last day, w1, w2, w3
days = [1 2 3];
roi_colors = max(0, brewermap(4,'Oranges'));
figure('Name', 'Post-training visual left ROI');
hold on;
set(gca,'ColorOrder',roi_colors)
AP_errorfill(timevec,[avg_last_visual_left_ROI avg_visual_left_ROI'],[sem_last_visual_left sem_visual_left'],roi_colors)
ylabel('Fluorescence')
xlabel('Time from stimulus onset')
legend({'', 'last training day', '', 'week 1', '', 'week 2', '', 'week 3'})
AM_prettyfig

% - frontal left ROIs 
% get day indexes for 113
animal_id = 3;
all_training_days = passive_master_activity(animal_id).day;
post_training_passive_days_mask = passive_master_activity(animal_id).post_training_passive_days_mask;
post_training_passive_day_idx = find(post_training_passive_days_mask);

frontal_left_avg_ROI = nan(length(post_training_passive_day_idx),length(timevec),length(animals));
for animal_id=3:length(animals)
    all_training_days = passive_master_activity(animal_id).day;
    post_training_passive_days_mask = passive_master_activity(animal_id).post_training_passive_days_mask;
    post_training_passive_day_idx = find(post_training_passive_days_mask);

    for day_idx=post_training_passive_day_idx
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
        deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
        
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        right_stim_act = deconvolved_stim_avg_act(:,:,trialStimulusValue==1);
        frontal_left = AP_svd_roi(U_master(:,:,1:num_comp),right_stim_act,[],[],roi.frontal_left);
        frontal_left = permute(frontal_left,[2,3,1]);
        frontal_left_avg_ROI(day_idx-post_training_passive_day_idx(1)+1,:,animal_id) = mean(frontal_left,2);
    end
end

% average across animals
this_frontal_left_avg_ROI = frontal_left_avg_ROI(:,:,3:end);
avg_frontal_left_ROI = nanmean(this_frontal_left_avg_ROI,3);
std_frontal_left = std(this_frontal_left_avg_ROI, 0, 3);
sem_frontal_left =  std_frontal_left./sqrt(length(animals)-2);


% get last day of training ROIs

last_day_training_idx = post_training_passive_day_idx(1)-1;

last_frontal_left_ROI = nan(length(timevec),length(animals));
day_idx = last_day_training_idx;
for animal_id=3:length(animals)
    stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
    deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
    deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
    
    trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
    right_stim_act = deconvolved_stim_avg_act(:,:,trialStimulusValue==1);
    frontal_left = AP_svd_roi(U_master(:,:,1:num_comp),right_stim_act,[],[],roi.frontal_left);
    frontal_left = permute(frontal_left,[2,3,1]);
    last_frontal_left_ROI(:,animal_id) = mean(frontal_left,2);
end

% average across animals
this_frontal_left_avg_ROI = last_frontal_left_ROI(:,3:end);
avg_last_frontal_left_ROI = nanmean(this_frontal_left_avg_ROI,2);
std_last_frontal_left = std(this_frontal_left_avg_ROI, 0, 2);
sem_last_frontal_left =  std_last_frontal_left./sqrt(length(animals)-2);

% make a plot of last day, w1, w2, w3
days = [1 2 3];
roi_colors = max(0, brewermap(4,'Reds')-0.2);
figure('Name', 'Post-training frontal left ROI');
hold on;
set(gca,'ColorOrder',roi_colors)
AP_errorfill(timevec,[avg_last_frontal_left_ROI avg_frontal_left_ROI'],[sem_last_frontal_left sem_frontal_left'],roi_colors)
ylabel('Fluorescence')
xlabel('Time from stimulus onset')
legend({'', 'last training day', '', 'week 1', '', 'week 2', '', 'week 3'})
AM_prettyfig

%% choose time window

% averaging window
small_window = [0.05 0.2];
small_window_idx = timevec>=small_window(1)&timevec<=small_window(2);

% frontal left
avg_frontal_left = nan(1,length(original_task_day_idx));
for day_idx=original_task_day_idx
    this_idx = day_idx-original_task_day_idx(1)+1;
    avg_frontal_left(this_idx) = mean(frontal_left_avg_ROI(this_idx,small_window_idx),2);
end

% frontal right
avg_frontal_right = nan(1,length(reversal_task_day_idx));
for day_idx=reversal_task_day_idx
    this_idx = day_idx-reversal_task_day_idx(1)+1;
    avg_frontal_right(this_idx) = mean(frontal_right_avg_ROI(this_idx,small_window_idx),2);
end

% plot average fluorescence
figure;
plot(1:length(original_task_day_idx), avg_frontal_left, '-o')
hold on;
plot(1:length(reversal_task_day_idx),avg_frontal_right, 'b-o')
title('AP107 mPFC avg passive fluorescence')
xlabel('Training day')
ylabel('Avg fluorescence 5-200 ms after stim onset')
legend({'Original', 'Reversal'})


% 
% % frontal left
% avg_frontal_left = nan(1,length(passive_master_activity.day));
% for day_idx=1:length(passive_master_activity.day)
%     stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
%     deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
%     deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
%     
%     trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
%     right_stim_act = deconvolved_stim_avg_act(:,:,trialStimulusValue==1);
%     frontal_left = AP_svd_roi(U_master(:,:,1:num_comp),right_stim_act,[],[],roi.frontal_left);
%     frontal_left = permute(frontal_left,[2,3,1]);
%     frontal_left_avg_ROI = mean(frontal_left(small_window,:),2);
%     avg_frontal_left(day_idx) = mean(frontal_left_avg_ROI,1);
% end
% 
% % frontal right
% avg_frontal_right = nan(1,length(new_task_days));
% for day_idx=1:length(new_task_days)
%     stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
%     deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
%     deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
%     trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
%     left_stim_act = deconvolved_stim_avg_act(:,:,trialStimulusValue==-1);
%     frontal_right = AP_svd_roi(U_master(:,:,1:num_comp),left_stim_act,[],[],roi.frontal_right);
%     frontal_right = permute(frontal_right,[2,3,1]);
%     frontal_right_avg_ROI = mean(frontal_right(small_window,:),2);
%     avg_frontal_right(day_idx) = mean(frontal_right_avg_ROI,1);
% end
% 
% % plot average fluorescence
% figure;
% plot(avg_frontal_left(4:end), '-o')
% hold on;
% plot(first_new_task_day:length(passive_master_activity.day)-3,avg_frontal_right, 'b-o')
% title('AP107 left mPFC avg fluorescence')
% xline(first_new_task_day, 'r')
% xlabel('Training day')
% ylabel('Avg fluorescence 200-500 ms after stim onset')

%% Avg picture with scroll

% load activity from passive
load('passive_master_activity.mat');
animals = {passive_master_activity.animal};

timestep = passive_master_activity.timestep;
start_time = passive_master_activity.start_time;
end_time = passive_master_activity.end_time;
stim_frame = passive_master_activity.stim_frame;
timevec = passive_master_activity.timevec;
num_comp = passive_master_activity.num_comp;

% load ROI masks
load('ROIs.mat');

% load master U
master_u_fn = 'D:\Andrada\Master project\widefield_alignment\U_master';
load(master_u_fn);

% choose averaging window
small_window = [0.05 0.3];
small_window_idx = timevec>=small_window(1)&timevec<=small_window(2);

% possible stims
possible_stimuli = [-1 0 1];

for animal_id=1:length(animals)
    passive_days = passive_master_activity(animal_id).day;
    muscimol_days = passive_master_activity(animal_id).muscimol_days;
    
    % initialize avg_act over window
    all_window_avg_act = nan(num_comp, length(passive_master_activity(animal_id).day), length(possible_stimuli));
    
    for day_idx=1:length(passive_days)
        day = passive_days(day_idx);
        
        % skip if it's a muscimol day
        is_muscimol = find(contains(muscimol_days,day));
        if is_muscimol
            continue
        end
        
        % get stim act for this animal and day
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        
        % skip if it's day with no imaging
        if isempty(stim_act)
            continue
        end
        
        % load trial information and wheel move
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        stim_wheel_move = passive_master_activity(animal_id).stim_wheel_move{day_idx};
        
        % find average activity for each stimulus
        all_stim_avg_act = nan(size(stim_act,1),size(stim_act,2),length(possible_stimuli));
        for stim_idx =1:length(possible_stimuli)
            no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
            this_stim_act = stim_act(:,:,no_move_trials&trialStimulusValue==possible_stimuli(stim_idx));
            all_stim_avg_act(:,:,stim_idx) = nanmean(this_stim_act,3);
        end
        all_stim_avg_act = all_stim_avg_act - all_stim_avg_act(:,stim_frame,:);
        
        % deconvolve
        deconvolved_all_stim_avg_act = AP_deconv_wf(all_stim_avg_act, [], 1/timestep);
        
        % average across chosen window
        window_avg_act = mean(deconvolved_all_stim_avg_act(:,small_window_idx,:),2);
        all_window_avg_act(:,day_idx,:) = window_avg_act;
    end
    
    % get fluoresence
    all_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),all_window_avg_act);
    passive_master_activity(animal_id).all_avg_fluorescence = all_avg_fluorescence;
    
    % video
    AP_image_scroll(all_avg_fluorescence,1:length(passive_master_activity(animal_id).day)); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)])
    axis image;
end

% save the avg_fluorescence pics in struct to load easily
% save('passive_master_activity.mat', 'passive_master_activity', '-v7.3')


%% Pictures each day + ROIs

figsdir = 'D:\Andrada\Master project\Random figures\';

% load activity from passive
load('passive_master_activity.mat');
animals = {passive_master_activity.animal};

timestep = passive_master_activity.timestep;
start_time = passive_master_activity.start_time;
end_time = passive_master_activity.end_time;
stim_frame = passive_master_activity.stim_frame;
timevec = passive_master_activity.timevec;
num_comp = passive_master_activity.num_comp;

% load ROI masks
load('ROIs.mat');

% load master U
master_u_fn = 'D:\Andrada\Master project\widefield_alignment\U_master';
load(master_u_fn);

% choose averaging window
small_window = [0.05 0.3];
small_window_idx = timevec>=small_window(1)&timevec<=small_window(2);

% possible stims
possible_stimuli = [-1 0 1];

for animal_id=1:length(animals)
    
    animal = animals(animal_id);
    
    % get fluoresence
    all_avg_fluorescence = passive_master_activity(animal_id).all_avg_fluorescence;
    
    avg_frontal_left = nan(3,length(passive_master_activity(animal_id).day));
    avg_frontal_right = nan(3,length(passive_master_activity(animal_id).day));
    
    for stim_idx =1:length(possible_stimuli)
        figure('Name',['Plot for animal ' cell2mat(animal) ' and stim ' num2str(possible_stimuli(stim_idx))]);
        for day_idx=1:length(passive_master_activity(animal_id).day)
            
            % load activity
            stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
            
            % skip if no activity
            if isempty(stim_act)
                continue
            end
            
            % brain pictures subplots
            subplot(5,10,day_idx);
            imagesc(all_avg_fluorescence(:,:,day_idx,stim_idx)); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)])
            title(['Day ' num2str(day_idx)])
            
            % deconvolve and baseline substract
            deconvolved_stim_avg_act = AP_deconv_wf(this_stim_act, [], 1/timestep);
            deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
            
            % get trial information
            trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
            
            % activity for current stim and no move
            no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
            this_stim_act = deconvolved_stim_avg_act(:,:,no_move_trials&trialStimulusValue==possible_stimuli(stim_idx));
            
            % frontal left ROI subplot
            frontal_left = AP_svd_roi(U_master(:,:,1:num_comp),this_stim_act,[],[],roi.frontal_left);
            frontal_left = permute(frontal_left,[2,3,1]);
            frontal_left_avg_ROI = mean(frontal_left(small_window_idx,:),2);
            avg_frontal_left(stim_idx,day_idx) = mean(frontal_left_avg_ROI,1);
            
            % frontal right ROI subplot
            frontal_right = AP_svd_roi(U_master(:,:,1:num_comp),this_stim_act,[],[],roi.frontal_right);
            frontal_right = permute(frontal_right,[2,3,1]);
            frontal_right_avg_ROI = mean(frontal_right(small_window_idx,:),2);
            avg_frontal_right(stim_idx,day_idx) = mean(frontal_right_avg_ROI,1);
        end
        
        % save in struct
        passive_master_activity(animal_id).avg_frontal_left = avg_frontal_left;
        passive_master_activity(animal_id).avg_frontal_right = avg_frontal_right;
        
        % subplots
        subplot(5,10,31:40);
        plot(avg_frontal_left(stim_idx,:), '-o');
        ylabel('Frontal left fluorescence');
        ylim([-1*10^(-3) 4*10^(-3)]);
        
        subplot(5,10,41:50);
        plot(avg_frontal_right(stim_idx,:), '-o');
        ylabel('Frontal right fluorescence');
        ylim([-1*10^(-3) 4*10^(-3)]);
        
        % looks bad because saves in small view
        % saveas(gcf,[figsdir 'Activity per day and ROIs for animal ' cell2mat(animal) ' and stim ' num2str(possible_stimuli(stim_idx)) '.png'])
    end
end

%% Different averaging windows

% load activity from passive
load('passive_master_activity.mat');
animals = {passive_master_activity.animal};

timestep = passive_master_activity.timestep;
start_time = passive_master_activity.start_time;
end_time = passive_master_activity.end_time;
stim_frame = passive_master_activity.stim_frame;
timevec = passive_master_activity.timevec;
num_comp = passive_master_activity.num_comp;

% load ROI masks
load('ROIs.mat');

% load master U
master_u_fn = 'D:\Andrada\Master project\widefield_alignment\U_master';
load(master_u_fn);

% choose averaging window
small_window = [0.05 0.2];
small_window_idx = timevec>=small_window(1)&timevec<=small_window(2);

% possible stims
possible_stimuli = [-1 0 1];

for animal_id=1:length(animals)
    passive_days = passive_master_activity(animal_id).day;
    muscimol_days = passive_master_activity(animal_id).muscimol_days;
    
    % initialize avg_act over window
    all_window_avg_act = nan(num_comp, length(passive_master_activity(animal_id).day), length(possible_stimuli));
    
    for day_idx=1:length(passive_days)
        day = passive_days(day_idx);
        
        % skip if it's a muscimol day
        is_muscimol = find(contains(muscimol_days,day));
        if is_muscimol
            continue
        end
        
        % get stim act for this animal and day
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        
        % skip if it's day with no imaging
        if isempty(stim_act)
            continue
        end
        
        % load trial information and wheel move
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        stim_wheel_move = passive_master_activity(animal_id).stim_wheel_move{day_idx};
        
        % find average activity for each stimulus
        all_stim_avg_act = nan(size(stim_act,1),size(stim_act,2),length(possible_stimuli));
        for stim_idx =1:length(possible_stimuli)
            no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
            this_stim_act = stim_act(:,:,no_move_trials&trialStimulusValue==possible_stimuli(stim_idx));
            
            %             % temp
            %             temp_deconvolved_act = AP_deconv_wf(this_stim_act, [], 1/timestep);
            %             temp_avg_deconvolved_act = mean(temp_deconvolved_act,3);
            %             temp_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),temp_avg_deconvolved_act);
            %             AP_image_scroll(temp_fluorescence,timevec); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)])
            %             axis image;
            %
            %             temp_stim_act = stim_act(:,:,trialStimulusValue==possible_stimuli(stim_idx));
            %             temp_deconvolved_act = AP_deconv_wf(temp_stim_act, [], 1/timestep);
            %             temp_avg_deconvolved_act = mean(temp_deconvolved_act,3);
            %             temp_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),temp_avg_deconvolved_act);
            %             AP_image_scroll(temp_fluorescence,timevec); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)])
            %             axis image;
            
            all_stim_avg_act(:,:,stim_idx) = nanmean(this_stim_act,3);
        end
        all_stim_avg_act = all_stim_avg_act - all_stim_avg_act(:,stim_frame,:);
        
        % deconvolve
        deconvolved_all_stim_avg_act = AP_deconv_wf(all_stim_avg_act, [], 1/timestep);
        
        % average across chosen window
        window_avg_act = mean(deconvolved_all_stim_avg_act(:,small_window_idx,:),2);
        all_window_avg_act(:,day_idx,:) = window_avg_act;
    end
    
    % get fluoresence
    all_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),all_window_avg_act);
    passive_master_activity(animal_id).all_avg_fluorescence_shorter = all_avg_fluorescence;
    
    % video
    AP_image_scroll(all_avg_fluorescence,1:length(passive_master_activity(animal_id).day)); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)])
    axis image;
end

%% ROI plots per task

figsdir = 'D:\Andrada\Master project\Random figures\Passive\ROI plots per task\';

% load activity from passive
load('passive_master_activity.mat');
animals = {passive_master_activity.animal};

timestep = passive_master_activity.timestep;
start_time = passive_master_activity.start_time;
end_time = passive_master_activity.end_time;
stim_frame = passive_master_activity.stim_frame;
timevec = passive_master_activity.timevec;
num_comp = passive_master_activity.num_comp;

% load ROI masks
load('ROIs.mat');
all_ROIs = fieldnames(roi);

% load master U
master_u_fn = 'D:\Andrada\Master project\widefield_alignment\U_master';
load(master_u_fn);

% choose averaging window
small_window = [0.05 0.2];
small_window_idx = timevec>=small_window(1)&timevec<=small_window(2);

% possible stims
possible_stimuli = [-1 0 1];

%% - Original task plot
% go through all ROIs
for ROI_idx=1:length(all_ROIs)
    
    this_ROI_values = roi.(all_ROIs{ROI_idx});
    
    % go through all stims
    for stim_idx =1:length(possible_stimuli)
        
        % name of plot
        figure('Name',['ROI trace ' cell2mat(all_ROIs(ROI_idx)) ' for stim ' num2str(possible_stimuli(stim_idx))]);
        
        % go through each animal
        for animal_id=1:length(animals)
            
            animal = animals(animal_id);
            
            % get day indexes
            all_training_days = passive_master_activity(animal_id).day;
            original_task_days_mask = passive_master_activity(animal_id).original_task_days_mask;
            muscimol_days_mask = passive_master_activity(animal_id).muscimol_days_mask;
            original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);
            first_day_idx = original_task_day_idx(1);
            
            % initialize avg_ROI trace to plot for this animal
            avg_this_ROI = nan(3,length(original_task_day_idx));
            
            for day_idx = original_task_day_idx
                
                % load activity
                stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
                
                % skip if no activity
                if isempty(stim_act)
                    continue
                end
                
                % deconvolve and baseline substract
                deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
                deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
                
                % get trial information
                trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
                stim_wheel_move = passive_master_activity(animal_id).stim_wheel_move{day_idx};
                
                % activity for current stim and no move
                no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
                this_stim_act = deconvolved_stim_avg_act(:,:,no_move_trials&trialStimulusValue==possible_stimuli(stim_idx));
                
                % get avg ROI
                this_ROI = AP_svd_roi(U_master(:,:,1:num_comp),this_stim_act,[],[],this_ROI_values);
                this_ROI = permute(this_ROI,[2,3,1]);
                this_avg_ROI = mean(this_ROI(small_window_idx,:),2);
                avg_this_ROI(stim_idx,day_idx-first_day_idx+1) = mean(this_avg_ROI,1);
            end
            
            % plot
            plot(avg_this_ROI(stim_idx,:), '-o', 'LineWidth', 2);
            hold on;
        end
        ylabel('ROI fluorescence');
        xlabel('Training day');
        ylim([-2*10^(-3) 6*10^(-3)]);
        legend(animals)
        AM_prettyfig
        % saveas(gcf,[figsdir 'Original_task_ROI trace ' cell2mat(all_ROIs(ROI_idx)) ' for stim ' num2str(possible_stimuli(stim_idx)) '.png'])
    end
end

%% - Reversal task plot

% make a colorOrder thing where 107's line is black and others are gray 

% colors_reversal = cat(1, [0 0 0], repmat([0.5 0.5 0.5], [4, 1]));

% go through all ROIs
for ROI_idx=1:length(all_ROIs)
    
    this_ROI_values = roi.(all_ROIs{ROI_idx});
    
    % go through all stims
    for stim_idx =1:length(possible_stimuli)
        
        % name of plot
        figure('Name',['ROI trace ' cell2mat(all_ROIs(ROI_idx)) ' for stim ' num2str(possible_stimuli(stim_idx))]);
        ax = axes;

        % go through each animal
        for animal_id=1:length(animals)
            
            animal = animals(animal_id);
            
            % get day indexes
            all_training_days = passive_master_activity(animal_id).day;
            reversal_task_days_mask = passive_master_activity(animal_id).reversal_task_days_mask;
            reversal_task_day_idx = find(reversal_task_days_mask);
            first_day_idx = reversal_task_day_idx(1);
            
            % initialize avg_ROI trace to plot for this animal
            avg_this_ROI = nan(3,length(reversal_task_day_idx));
            
            for day_idx = reversal_task_day_idx
                
                % load activity
                stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
                
                % skip if no activity
                if isempty(stim_act)
                    continue
                end
                
                % deconvolve and baseline substract
                deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
                deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
                
                % get trial information
                trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
                stim_wheel_move = passive_master_activity(animal_id).stim_wheel_move{day_idx};
                
                % activity for current stim and no move
                no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
                this_stim_act = deconvolved_stim_avg_act(:,:,no_move_trials&trialStimulusValue==possible_stimuli(stim_idx));
                
                % get avg ROI
                this_ROI = AP_svd_roi(U_master(:,:,1:num_comp),this_stim_act,[],[],this_ROI_values);
                this_ROI = permute(this_ROI,[2,3,1]);
                this_avg_ROI = mean(this_ROI(small_window_idx,:),2);
                avg_this_ROI(stim_idx,day_idx-first_day_idx+1) = mean(this_avg_ROI,1);
            end
            
            % plot
            plot(avg_this_ROI(stim_idx,:), '-o', 'LineWidth', 2);
            hold on;
        end
        ylabel('ROI fluorescence');
        xlabel('Training day');
        ylim([-2*10^(-3) 6*10^(-3)]);
%         ax.ColorOrder = colors_reversal;
%         legend({'learner', 'non-learners'})
        legend(animals)
        AM_prettyfig
        % saveas(gcf,[figsdir 'Reversal_task_ROI trace ' cell2mat(all_ROIs(ROI_idx)) ' for stim ' num2str(possible_stimuli(stim_idx)) '.png'])
    end
end

%% Average brain pictures

figsdir = 'D:\Andrada\Master project\Random figures\ROI plots per task\';

% load activity from passive
load('passive_master_activity.mat');
animals = {passive_master_activity.animal};

timestep = passive_master_activity.timestep;
start_time = passive_master_activity.start_time;
end_time = passive_master_activity.end_time;
stim_frame = passive_master_activity.stim_frame;
timevec = passive_master_activity.timevec;
num_comp = passive_master_activity.num_comp;

% load ROI masks
load('ROIs.mat');
all_ROIs = fieldnames(roi);

% load master U
master_u_fn = 'D:\Andrada\Master project\widefield_alignment\U_master';
load(master_u_fn);

% choose averaging window
small_window = [0.05 0.2];
small_window_idx = timevec>=small_window(1)&timevec<=small_window(2);

% possible stims
possible_stimuli = [-1 0 1];

%% - Passive

% temp averaging window
% small_window = [0.05 0.3];
% small_window_idx = timevec>=small_window(1)&timevec<=small_window(2);

% learned day in task - passive day 1 (should all days but didn't have time to think)
all_passive_days = 1:3;

% to have dim for avg activity from size of first master activity
stim_act = passive_master_activity(1).stim_activity{1};
all_stim_avg_act = nan(size(stim_act,1),size(stim_act,2),length(possible_stimuli));

for stim_idx =1:length(possible_stimuli)
    this_stim_act = [];
    for passive_day = all_passive_days
        for animal_id=1:length(animals)
            
            % get day indexes
            all_training_days = passive_master_activity(animal_id).day;
            passive_only_days_mask = passive_master_activity(animal_id).passive_only_days_mask;
            passive_only_day_idx = find(passive_only_days_mask);
            first_day_idx = passive_only_day_idx(1);
            
            % find which day the learned day is in whole dataset
            day_idx = passive_only_day_idx(passive_day);
            day = all_training_days{day_idx};
            
            % load and avg
            stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
            trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
            stim_wheel_move = passive_master_activity(animal_id).stim_wheel_move{day_idx};
            no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
            this_mouse_stim_act = stim_act(:,:,no_move_trials&trialStimulusValue==possible_stimuli(stim_idx));
            this_stim_act = cat(3, this_stim_act, this_mouse_stim_act);
        end
    end
    all_stim_avg_act(:,:,stim_idx) = nanmean(this_stim_act,3);
end
all_stim_avg_act = all_stim_avg_act - all_stim_avg_act(:,stim_frame,:);

deconvolved_all_stim_avg_act = AP_deconv_wf(all_stim_avg_act, [], 1/timestep);

% also average over window
passive_only_avg_window_act = mean(deconvolved_all_stim_avg_act(:,small_window_idx,:),2);
passive_only_avg_window_act = permute(passive_only_avg_window_act, [1, 3, 2]);

% get fluoresence
all_stim_interval_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),passive_only_avg_window_act);

% make subplot for each stim
figure; title('Passive all days');
for stim_idx =1:length(possible_stimuli)
    subplot(3,1,stim_idx);
    imagesc(all_stim_interval_avg_fluorescence(:,:,stim_idx)); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)]);
    axis image;
    set(gca,'Xtick',[])
    set(gca,'Ytick',[])
end

% plot for each stim
for stim_idx =1:length(possible_stimuli)
    figure; title(['Passive all days' num2str(possible_stimuli(stim_idx))]);
    imagesc(all_stim_interval_avg_fluorescence(:,:,stim_idx)); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)]);
    axis image;
    axis off;
    hold on;
    AP_reference_outline('ccf_aligned','#808080');
end
%% - Original learned day

% learned day in task
learned_day_idx = 5;

% to have dim for avg activity from size of first master activity
stim_act = passive_master_activity(1).stim_activity{1};
all_stim_avg_act = nan(size(stim_act,1),size(stim_act,2),length(possible_stimuli));

for stim_idx =1:length(possible_stimuli)
    this_stim_act = [];
    for animal_id=1:length(animals)
        
        % get day indexes
        all_training_days = passive_master_activity(animal_id).day;
        original_task_days_mask = passive_master_activity(animal_id).original_task_days_mask;
        original_task_day_idx = find(original_task_days_mask);
        first_day_idx = original_task_day_idx(1);
        
        % find which day the learned day is in whole dataset
        these_day_idx = original_task_day_idx(learned_day_idx:end);
        
        % go through all learned days
        for day_idx = these_day_idx
            day = all_training_days{day_idx};
            
            % load and avg
            stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
            trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
            stim_wheel_move = passive_master_activity(animal_id).stim_wheel_move{day_idx};
            no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
            this_mouse_stim_act = stim_act(:,:,no_move_trials&trialStimulusValue==possible_stimuli(stim_idx));
            this_stim_act = cat(3, this_stim_act, this_mouse_stim_act);
        end
    end
    all_stim_avg_act(:,:,stim_idx) = nanmean(this_stim_act,3);
end
all_stim_avg_act = all_stim_avg_act - all_stim_avg_act(:,stim_frame,:);

deconvolved_all_stim_avg_act = AP_deconv_wf(all_stim_avg_act, [], 1/timestep);

% also average over window
original_task_avg_window_act = mean(deconvolved_all_stim_avg_act(:,small_window_idx,:),2);
original_task_avg_window_act = permute(original_task_avg_window_act, [1, 3, 2]);

% get fluoresence
all_stim_interval_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),original_task_avg_window_act);

% make subplot for each stim
figure; title(['Original all from ' num2str(learned_day_idx)]);
for stim_idx =1:length(possible_stimuli)
    subplot(3,1,stim_idx);
    imagesc(all_stim_interval_avg_fluorescence(:,:,stim_idx)); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)]);
    axis image;
    set(gca,'Xtick',[])
    set(gca,'Ytick',[])
end

% plot for each stim
for stim_idx =1:length(possible_stimuli)
    figure('Name',['Original task all post-learning days stim' num2str(possible_stimuli(stim_idx))]);
    imagesc(all_stim_interval_avg_fluorescence(:,:,stim_idx)); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)]);
    axis image;
    axis off;
    hold on;
    AP_reference_outline('ccf_aligned','#808080');
end
%% - Reversal learned day

% learned day in task
learned_day_idx = 7;

% to have dim for avg activity from size of first master activity
stim_act = passive_master_activity(1).stim_activity{1};
all_stim_avg_act = nan(size(stim_act,1),size(stim_act,2),length(possible_stimuli));

%% -- AP107 - learner
animal_id = 1;

% get day indexes
all_training_days = passive_master_activity(animal_id).day;
reversal_task_days_mask = passive_master_activity(animal_id).reversal_task_days_mask;
reversal_task_day_idx = find(reversal_task_days_mask);
first_day_idx = reversal_task_day_idx(1);

% find day indices for learned days
these_day_idx = reversal_task_day_idx(learned_day_idx:end);

% go through all learned days
for stim_idx =1:length(possible_stimuli)
    this_stim_act = [];
    
    for day_idx = these_day_idx
        day = all_training_days{day_idx};
        
        % load and avg
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        stim_wheel_move = passive_master_activity(animal_id).stim_wheel_move{day_idx};
        no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
        this_mouse_stim_act = stim_act(:,:,no_move_trials&trialStimulusValue==possible_stimuli(stim_idx));
        this_stim_act = cat(3, this_stim_act, this_mouse_stim_act);
        
        % average for stim
        all_stim_avg_act(:,:,stim_idx) = nanmean(this_stim_act,3);
    end
    
end
% substract baseline
all_stim_avg_act = all_stim_avg_act - all_stim_avg_act(:,stim_frame,:);

% deconvolve
deconvolved_all_stim_avg_act = AP_deconv_wf(all_stim_avg_act, [], 1/timestep);

% also average over window
reversal_task_avg_window_act = mean(deconvolved_all_stim_avg_act(:,small_window_idx,:),2);
reversal_task_avg_window_act = permute(reversal_task_avg_window_act, [1, 3, 2]);

% get fluoresence
all_stim_interval_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),reversal_task_avg_window_act);

% make subplot for each stim
figure; title(['Reversal learner from ' num2str(learned_day_idx)]);
for stim_idx =1:length(possible_stimuli)
    subplot(3,1,stim_idx);
    imagesc(all_stim_interval_avg_fluorescence(:,:,stim_idx)); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)]);
    axis image;
    set(gca,'Xtick',[])
    set(gca,'Ytick',[])
end

for stim_idx =1:length(possible_stimuli)
    figure('Name',['Reversal task all post-learning days for learner stim' num2str(possible_stimuli(stim_idx))]);
    imagesc(all_stim_interval_avg_fluorescence(:,:,stim_idx)); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)]);
    axis image;
    axis off;
    hold on;
    AP_reference_outline('ccf_aligned','#808080');
end

%% -- other mice

for stim_idx =1:length(possible_stimuli)
    this_stim_act = [];
    for animal_id=1:length(animals)
        
        % skip 107
        if animal_id == 1
            continue
        end
        
        % get day indexes
        all_training_days = passive_master_activity(animal_id).day;
        reversal_task_days_mask = passive_master_activity(animal_id).reversal_task_days_mask;
        reversal_task_day_idx = find(reversal_task_days_mask);
        first_day_idx = reversal_task_day_idx(1);
        
        % find day indices for learned days
        these_day_idx = reversal_task_day_idx(learned_day_idx:end);
        
        for day_idx = these_day_idx
            day = all_training_days{day_idx};
            
            % load and avg
            stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
            trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
            stim_wheel_move = passive_master_activity(animal_id).stim_wheel_move{day_idx};
            no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
            this_mouse_stim_act = stim_act(:,:,no_move_trials&trialStimulusValue==possible_stimuli(stim_idx));
            this_stim_act = cat(3, this_stim_act, this_mouse_stim_act);
        end
        all_stim_avg_act(:,:,stim_idx) = nanmean(this_stim_act,3);
    end
end

all_stim_avg_act = all_stim_avg_act - all_stim_avg_act(:,stim_frame,:);

deconvolved_all_stim_avg_act = AP_deconv_wf(all_stim_avg_act, [], 1/timestep);

% also average over window
reversal_task_avg_window_act = mean(deconvolved_all_stim_avg_act(:,small_window_idx,:),2);
reversal_task_avg_window_act = permute(reversal_task_avg_window_act, [1, 3, 2]);

% get fluoresence
all_stim_interval_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),reversal_task_avg_window_act);

% make subplot for each stim
figure; title(['Reversal non-learners from ' num2str(learned_day_idx)]);
for stim_idx =1:length(possible_stimuli)
    subplot(3,1,stim_idx);
    imagesc(all_stim_interval_avg_fluorescence(:,:,stim_idx)); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)]);
    axis image;
    set(gca,'Xtick',[])
    set(gca,'Ytick',[])
end

for stim_idx =1:length(possible_stimuli)
    figure('Name',['Reversal task all post-learning days for non-learners stim' num2str(possible_stimuli(stim_idx))]);
    imagesc(all_stim_interval_avg_fluorescence(:,:,stim_idx)); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)]);
    axis image;
    axis off;
    hold on;
    AP_reference_outline('ccf_aligned','#808080');
end


%% - Post-training
% load activity from passive
load('passive_master_activity.mat');
animals = {passive_master_activity.animal};

timestep = passive_master_activity.timestep;
start_time = passive_master_activity.start_time;
end_time = passive_master_activity.end_time;
stim_frame = passive_master_activity.stim_frame;
timevec = passive_master_activity.timevec;
num_comp = passive_master_activity.num_comp;

% load ROI masks
load('ROIs.mat');
all_ROIs = fieldnames(roi);

% load master U
master_u_fn = 'D:\Andrada\Master project\widefield_alignment\U_master';
load(master_u_fn);

% choose averaging window
small_window = [0.05 0.2];
small_window_idx = timevec>=small_window(1)&timevec<=small_window(2);

% possible stims
possible_stimuli = [-1 0 1];

% to have dim for avg activity from size of first master activity
stim_act = passive_master_activity(1).stim_activity{1};
all_stim_avg_act = nan(size(stim_act,1),size(stim_act,2),length(possible_stimuli), 3);

% get day indexes for 115
animal_id = 5;
all_training_days = passive_master_activity(animal_id).day;
post_training_passive_days_mask = passive_master_activity(animal_id).post_training_passive_days_mask;
post_training_passive_days_idx = find(post_training_passive_days_mask);

% day indices
these_day_idx = post_training_passive_days_idx;

for stim_idx =1:length(possible_stimuli)
    for day_idx = these_day_idx
        this_stim_act = [];
        for animal_id=3:length(animals) % skip 107 and 108
   
            % load and avg
            stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
            trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
            stim_wheel_move = passive_master_activity(animal_id).stim_wheel_move{day_idx};
            no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
            this_mouse_stim_act = stim_act(:,:,no_move_trials&trialStimulusValue==possible_stimuli(stim_idx));
            this_stim_act = cat(3, this_stim_act, this_mouse_stim_act);
        end
        all_stim_avg_act(:,:,stim_idx,day_idx-these_day_idx(1)+1) = nanmean(this_stim_act,3);
    end
end

all_stim_avg_act = all_stim_avg_act - all_stim_avg_act(:,stim_frame,:,:);

deconvolved_all_stim_avg_act = AP_deconv_wf(all_stim_avg_act, [], 1/timestep);

% also average over window
passive_only_avg_window_act = mean(deconvolved_all_stim_avg_act(:,small_window_idx,:,:),2);
passive_only_avg_window_act = permute(passive_only_avg_window_act, [1, 3, 4, 2]);

% plot for each stim
for day_idx = 1:3
    
    this_act = passive_only_avg_window_act(:,:,day_idx);
    % get fluoresence
    all_stim_interval_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),this_act);
    
    for stim_idx =1:length(possible_stimuli)
        figure('Name', ['Passive day ' num2str(day_idx) ' stim' num2str(possible_stimuli(stim_idx))]);
        imagesc(all_stim_interval_avg_fluorescence(:,:,stim_idx)); colormap(brewermap([], 'PRGn'));caxis([-3.5*10^(-3) 3.5*10^(-3)]);
        axis image;
        axis off;
        hold on;
        AP_reference_outline('ccf_aligned','#808080');
    end
end

%% -- get last day of training

% to have dim for avg activity from size of first master activity
stim_act = passive_master_activity(1).stim_activity{1};
all_stim_avg_act = nan(size(stim_act,1),size(stim_act,2),length(possible_stimuli));

% get day
animal_id = 5;
all_training_days = passive_master_activity(animal_id).day;
post_training_passive_days_mask = passive_master_activity(animal_id).post_training_passive_days_mask;
post_training_passive_days_idx = find(post_training_passive_days_mask);

% day indices
these_day_idx = post_training_passive_days_idx-1;

day_idx = post_training_passive_days_idx(1)-1;

for stim_idx =1:length(possible_stimuli)
    this_stim_act = [];
    for animal_id=3:length(animals) % skip 107 and 108
        
        % load and avg
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        stim_wheel_move = passive_master_activity(animal_id).stim_wheel_move{day_idx};
        no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
        this_mouse_stim_act = stim_act(:,:,no_move_trials&trialStimulusValue==possible_stimuli(stim_idx));
        this_stim_act = cat(3, this_stim_act, this_mouse_stim_act);
    end
    all_stim_avg_act(:,:,stim_idx) = nanmean(this_stim_act,3);
end

all_stim_avg_act = all_stim_avg_act - all_stim_avg_act(:,stim_frame,:,:);

deconvolved_all_stim_avg_act = AP_deconv_wf(all_stim_avg_act, [], 1/timestep);

% also average over window
passive_only_avg_window_act = mean(deconvolved_all_stim_avg_act(:,small_window_idx,:),2);
passive_only_avg_window_act = permute(passive_only_avg_window_act, [1, 3, 2]);

% get fluoresence
all_stim_interval_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),passive_only_avg_window_act);

for stim_idx =1:length(possible_stimuli)
    figure('Name', ['Passive day ' num2str(day_idx) ' stim' num2str(possible_stimuli(stim_idx))]);
    imagesc(all_stim_interval_avg_fluorescence(:,:,stim_idx)); colormap(brewermap([], 'PRGn'));caxis([-3.5*10^(-3) 3.5*10^(-3)]);
    axis image;
    axis off;
    hold on;
    AP_reference_outline('ccf_aligned','#808080');
end


%% Stim and ROI plot

figsdir = 'D:\Andrada\Master project\Random figures\ROI plots per task\';

% load activity from passive
load('passive_master_activity.mat');
animals = {passive_master_activity.animal};

timestep = passive_master_activity.timestep;
start_time = passive_master_activity.start_time;
end_time = passive_master_activity.end_time;
stim_frame = passive_master_activity.stim_frame;
timevec = passive_master_activity.timevec;
num_comp = passive_master_activity.num_comp;

% load ROI masks
load('ROIs.mat');
all_ROIs = fieldnames(roi);

% load master U
master_u_fn = 'D:\Andrada\Master project\widefield_alignment\U_master';
load(master_u_fn);

% choose averaging window
small_window = [0.05 0.2];
small_window_idx = timevec>=small_window(1)&timevec<=small_window(2);

% possible stims
possible_stimuli = [-1 0 1];

%% - original task

% define learned day
learned_day_idx = 5;

% initialize frontal ROIs
all_left_frontal_ROIs = nan(length(animals),length(timevec));
all_right_frontal_ROIs = nan(length(animals),length(timevec));

for animal_id=1:length(animals)
    
    % get day indexes
    all_training_days = passive_master_activity(animal_id).day;
    original_task_days_mask = passive_master_activity(animal_id).original_task_days_mask;
    original_task_day_idx = find(original_task_days_mask);
    
    % find which day the learned day is in whole dataset
    these_day_idx = original_task_day_idx(learned_day_idx:end);
    
    % initialize empty ROIs thing
    this_left_ROI = nan(length(timevec),length(these_day_idx));
    this_right_ROI = nan(length(timevec),length(these_day_idx));
    
    % go through all learned days
    for day_idx = these_day_idx
        day = all_training_days{day_idx};
        
        % load activity
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        
        % deconvolve and baseline substract
        deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
        deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
        
        % get trial information
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        stim_wheel_move = passive_master_activity(animal_id).stim_wheel_move{day_idx};
        
        % trials with no move
        no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
        
        % right stim act
        right_stim_act = deconvolved_stim_avg_act(:,:,no_move_trials&trialStimulusValue==1);
        left_stim_act = deconvolved_stim_avg_act(:,:,no_move_trials&trialStimulusValue==-1);
        
        % get avg left ROI for right stim
        left_ROI = AP_svd_roi(U_master(:,:,1:num_comp),right_stim_act,[],[],roi.frontal_left);
        left_ROI = permute(left_ROI,[2,3,1]);
        this_left_ROI(:,day_idx-these_day_idx(1)+1) = mean(left_ROI,2);
        
        % get avg right ROI for left stim
        right_ROI = AP_svd_roi(U_master(:,:,1:num_comp),left_stim_act,[],[],roi.frontal_right);
        right_ROI = permute(right_ROI,[2,3,1]);
        this_right_ROI(:,day_idx-these_day_idx(1)+1) = mean(right_ROI,2);
    end
    
    all_left_frontal_ROIs(animal_id,:) = mean(this_left_ROI,2);
    all_right_frontal_ROIs(animal_id,:) = mean(this_right_ROI,2);
end

% plot
roi_colors = cat(1,brewermap(length(animals),'Reds'),brewermap(length(animals),'Blues'));
figure; title('Original task');
hold on;
set(gca,'ColorOrder',roi_colors)
plot(timevec, all_left_frontal_ROIs')
hold on;
plot(timevec, all_right_frontal_ROIs');
legend(animals)

% avg plot

avg_all_left_frontal_ROIs = mean(all_left_frontal_ROIs, 1);
avg_all_right_frontal_ROIs = mean(all_right_frontal_ROIs, 1);

roi_colors = cat(1,brewermap(1,'Reds'),brewermap(1,'Blues'));
figure; title('Original task all average');
hold on;
set(gca,'ColorOrder',roi_colors)
plot(timevec, avg_all_left_frontal_ROIs')
hold on;
plot(timevec, avg_all_right_frontal_ROIs');
xline(0)
ylim([-1.5*10^-3 3*10^-3])
legend({'Left','Right'})

% avg plot for non-learners

avg_all_left_frontal_ROIs = mean(all_left_frontal_ROIs(2:end,:), 1);
avg_all_right_frontal_ROIs = mean(all_right_frontal_ROIs(2:end,:), 1);

roi_colors = cat(1,brewermap(1,'Reds'),brewermap(1,'Blues'));
figure; title('Original task reversal non-learners');
hold on;
set(gca,'ColorOrder',roi_colors)
plot(timevec, avg_all_left_frontal_ROIs')
hold on;
plot(timevec, avg_all_right_frontal_ROIs');
xline(0)
ylim([-1.5*10^-3 3*10^-3])
legend({'Left','Right'})

% plot for learner
roi_colors = cat(1,brewermap(1,'Reds'),brewermap(1,'Blues'));
figure; title('Original task reversal learner');
hold on;
set(gca,'ColorOrder',roi_colors)
plot(timevec, all_left_frontal_ROIs(1,:)')
hold on;
plot(timevec, all_right_frontal_ROIs(1,:)');
xline(0)
ylim([-1.5*10^-3 3*10^-3])
legend({'Left','Right'})

% plot both on same plot
roi_colors = cat(1,brewermap(2,'Reds'),brewermap(2,'Blues'));
figure; title('Original task reversal learner and non-learners');
hold on;
set(gca,'ColorOrder',roi_colors)
plot(timevec, all_left_frontal_ROIs(1,:)')
hold on;
plot(timevec, avg_all_left_frontal_ROIs')
hold on;
plot(timevec, all_right_frontal_ROIs(1,:)');
hold on;
plot(timevec, avg_all_right_frontal_ROIs');
xline(0)
ylim([-1.5*10^-3 3*10^-3])
legend({'Left learner','Left non-learners', 'Right learner', 'Right non-learners'})


%% - reversal

% define learned day
learned_day_idx = 7;

% initialize frontal ROIs
all_left_frontal_ROIs = nan(length(animals),length(timevec));
all_right_frontal_ROIs = nan(length(animals),length(timevec));

for animal_id=1:length(animals)
    
    % get day indexes
    all_training_days = passive_master_activity(animal_id).day;
    reversal_task_days_mask = passive_master_activity(animal_id).reversal_task_days_mask;
    reversal_task_day_idx = find(reversal_task_days_mask);
    first_day_idx = reversal_task_day_idx(1);
    
    % find day indices for learned days
    these_day_idx = reversal_task_day_idx(learned_day_idx:end);
    
    % initialize empty ROIs thing
    this_left_ROI = nan(length(timevec),length(these_day_idx));
    this_right_ROI = nan(length(timevec),length(these_day_idx));
    
    % go through all learned days
    for day_idx = these_day_idx
        day = all_training_days{day_idx};
        
        % load activity
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        
        % deconvolve and baseline substract
        deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
        deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
        
        % get trial information
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        stim_wheel_move = passive_master_activity(animal_id).stim_wheel_move{day_idx};
        
        % trials with no move
        no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
        
        % right stim act
        right_stim_act = deconvolved_stim_avg_act(:,:,no_move_trials&trialStimulusValue==1);
        left_stim_act = deconvolved_stim_avg_act(:,:,no_move_trials&trialStimulusValue==-1);
        
        % get avg left ROI for right stim
        left_ROI = AP_svd_roi(U_master(:,:,1:num_comp),right_stim_act,[],[],roi.frontal_left);
        left_ROI = permute(left_ROI,[2,3,1]);
        this_left_ROI(:,day_idx-these_day_idx(1)+1) = mean(left_ROI,2);
        
        % get avg right ROI for left stim
        right_ROI = AP_svd_roi(U_master(:,:,1:num_comp),left_stim_act,[],[],roi.frontal_right);
        right_ROI = permute(right_ROI,[2,3,1]);
        this_right_ROI(:,day_idx-these_day_idx(1)+1) = mean(right_ROI,2);
    end
    
    all_left_frontal_ROIs(animal_id,:) = mean(this_left_ROI,2);
    all_right_frontal_ROIs(animal_id,:) = mean(this_right_ROI,2);
end

% plot
roi_colors = cat(1,brewermap(length(animals),'Reds'),brewermap(length(animals),'Blues'));
figure; title('Reversal task');
hold on;
set(gca,'ColorOrder',roi_colors)
plot(timevec, all_left_frontal_ROIs')
hold on;
plot(timevec, all_right_frontal_ROIs');
legend(animals)

% avg plot for non-learners

avg_all_left_frontal_ROIs = mean(all_left_frontal_ROIs(2:end,:), 1);
avg_all_right_frontal_ROIs = mean(all_right_frontal_ROIs(2:end,:), 1);

roi_colors = cat(1,brewermap(1,'Reds'),brewermap(1,'Blues'));
figure; title('Reversal task non-learners');
hold on;
set(gca,'ColorOrder',roi_colors)
plot(timevec, avg_all_left_frontal_ROIs')
hold on;
plot(timevec, avg_all_right_frontal_ROIs');
xline(0)
ylim([-1.5*10^-3 3*10^-3])
legend({'Left','Right'})

% plot for learner
roi_colors = cat(1,brewermap(1,'Reds'),brewermap(1,'Blues'));
figure; title('Reversal task learner');
hold on;
set(gca,'ColorOrder',roi_colors)
plot(timevec, all_left_frontal_ROIs(1,:)')
hold on;
plot(timevec, all_right_frontal_ROIs(1,:)');
xline(0)
ylim([-1.5*10^-3 3*10^-3])
legend({'Left','Right'})

% plot both on same plot
roi_colors = cat(1,brewermap(2,'Reds'),brewermap(2,'Blues'));
figure; title('Reversal task learner and non-learners');
hold on;
set(gca,'ColorOrder',roi_colors)
plot(timevec, all_left_frontal_ROIs(1,:)')
hold on;
plot(timevec, avg_all_left_frontal_ROIs')
hold on;
plot(timevec, all_right_frontal_ROIs(1,:)');
hold on;
plot(timevec, avg_all_right_frontal_ROIs');
xline(0)
ylim([-1.5*10^-3 3*10^-3])
legend({'Left learner','Left non-learners', 'Right learner', 'Right non-learners'})


%% - AP107 like task ROI plot

% load activity from passive
load('passive_master_activity.mat');
animals = {passive_master_activity.animal};

timestep = passive_master_activity.timestep;
start_time = passive_master_activity.start_time;
end_time = passive_master_activity.end_time;
stim_frame = passive_master_activity.stim_frame;
timevec = passive_master_activity.timevec;
num_comp = passive_master_activity.num_comp;

% load ROI masks
load('ROIs.mat');
all_ROIs = fieldnames(roi);

% load master U
master_u_fn = 'D:\Andrada\Master project\widefield_alignment\U_master';
load(master_u_fn);

% choose averaging window
small_window = [0.05 0.2];
small_window_idx = timevec>=small_window(1)&timevec<=small_window(2);

% possible stims
possible_stimuli = [-1 0 1];

%% original task

% define learned day
learned_day_idx = 5;

% initialize frontal ROIs
original_left_frontal_ROIs = nan(length(animals),length(timevec));
original_right_frontal_ROIs = nan(length(animals),length(timevec));

for animal_id=1:length(animals)
    
    % get day indexes
    all_training_days = passive_master_activity(animal_id).day;
    original_task_days_mask = passive_master_activity(animal_id).original_task_days_mask;
    muscimol_days_mask = passive_master_activity(animal_id).muscimol_days_mask;
    original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);
    
    % find which day the learned day is in whole dataset
    these_day_idx = original_task_day_idx(learned_day_idx:end);
    
    % initialize empty ROIs thing
    this_left_ROI = nan(length(timevec),length(these_day_idx));
    this_right_ROI = nan(length(timevec),length(these_day_idx));
    
    % go through all learned days
    for day_idx = these_day_idx
        day = all_training_days{day_idx};
        
        % load activity
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        
        % deconvolve and baseline substract
        deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
        deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
        
        % take no move trials and identify stim_id = 1
        % get trial information
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        stim_wheel_move = passive_master_activity(animal_id).stim_wheel_move{day_idx};
        
        % trials with no move
        no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
        
        % right stim act
        right_stim_act = deconvolved_stim_avg_act(:,:,no_move_trials&trialStimulusValue==1);
        
        % get avg left ROI for right stim
        left_ROI = AP_svd_roi(U_master(:,:,1:num_comp),right_stim_act,[],[],roi.frontal_left);
        left_ROI = permute(left_ROI,[2,3,1]);
        this_left_ROI(:,day_idx-these_day_idx(1)+1) = mean(left_ROI,2);
        
        % get avg right ROI for right stim
        right_ROI = AP_svd_roi(U_master(:,:,1:num_comp),right_stim_act,[],[],roi.frontal_right);
        right_ROI = permute(right_ROI,[2,3,1]);
        this_right_ROI(:,day_idx-these_day_idx(1)+1) = mean(right_ROI,2);
    end
    
    original_left_frontal_ROIs(animal_id,:) = mean(this_left_ROI,2);
    original_right_frontal_ROIs(animal_id,:) = mean(this_right_ROI,2);
end

% reversal task

% define learned day
learned_day_idx = 7;

% initialize frontal ROIs
reversal_left_frontal_ROIs = nan(length(animals),length(timevec));
reversal_right_frontal_ROIs = nan(length(animals),length(timevec));

for animal_id=1:length(animals)
    
    % get day indexes
    all_training_days = passive_master_activity(animal_id).day;
    reversal_task_days_mask = passive_master_activity(animal_id).reversal_task_days_mask;
    reversal_task_day_idx = find(reversal_task_days_mask);
    
    % find which day the learned day is in whole dataset
    these_day_idx = reversal_task_day_idx(learned_day_idx:end);
    
    % initialize empty ROIs thing
    this_left_ROI = nan(length(timevec),length(these_day_idx));
    this_right_ROI = nan(length(timevec),length(these_day_idx));
    
    % go through all learned days
    for day_idx = these_day_idx
        day = all_training_days{day_idx};
        
        % load activity
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        
        % deconvolve and baseline substract
        deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
        deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
        
        % take no move trials and identify stim_id = -1
        % get trial information
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        stim_wheel_move = passive_master_activity(animal_id).stim_wheel_move{day_idx};
        
        % trials with no move
        no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
        
        % left stim act
        left_stim_act = deconvolved_stim_avg_act(:,:,no_move_trials&trialStimulusValue==-1);
        
        % get avg left ROI for left stim
        left_ROI = AP_svd_roi(U_master(:,:,1:num_comp),left_stim_act,[],[],roi.frontal_left);
        left_ROI = permute(left_ROI,[2,3,1]);
        this_left_ROI(:,day_idx-these_day_idx(1)+1) = mean(left_ROI,2);
        
        % get avg right ROI for left stim
        right_ROI = AP_svd_roi(U_master(:,:,1:num_comp),left_stim_act,[],[],roi.frontal_right);
        right_ROI = permute(right_ROI,[2,3,1]);
        this_right_ROI(:,day_idx-these_day_idx(1)+1) = mean(right_ROI,2);
    end
    
    reversal_left_frontal_ROIs(animal_id,:) = mean(this_left_ROI,2);
    reversal_right_frontal_ROIs(animal_id,:) = mean(this_right_ROI,2);
end


% plot
roi_colors = cat(1,brewermap(2,'Reds'),brewermap(2,'Blues'));

figure;
title('mPFC ROIs in original and reversal passive learner')
hold on;
set(gca,'ColorOrder',roi_colors)
plot(timevec, original_left_frontal_ROIs(1,:)')
hold on;
plot(timevec, reversal_left_frontal_ROIs(1,:)');
hold on;
plot(timevec, original_right_frontal_ROIs(1,:)')
hold on;
plot(timevec, reversal_right_frontal_ROIs(1,:)');
hold on;
xline(0.1, '-', {'RT original'}, 'LineWidth', 1.5)
xline(0.2, '-', {'RT reversal'}, 'LineWidth', 3)
legend({'Original task left mPFC', 'Reversal task left mPFC', 'Original task right mPFC', 'Reversal task right mPFC'})



%% Averaged vis+frontal avg fluorescence 

% averaging window
small_window = [0.05 0.2];
small_window_idx = timevec>=small_window(1)&timevec<=small_window(2);

possible_stimuli = [-1 0 1];
%% - Original 

%% -- visual left
% stim left
stim_idx = 3;

% just initialize with nans for all animals and days of 107
% get day indexes
animal_id=1;
all_training_days = passive_master_activity(animal_id).day;
original_task_days_mask = passive_master_activity(animal_id).original_task_days_mask;
muscimol_days_mask = passive_master_activity(animal_id).muscimol_days_mask;
original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);
num_original_days = length(original_task_day_idx);

all_avg_visual_left = nan(length(animals), length(original_task_day_idx));

% go through each animal
for animal_id=1:length(animals)
    
    animal = animals(animal_id);
    
    % get day indexes
    all_training_days = passive_master_activity(animal_id).day;
    original_task_days_mask = passive_master_activity(animal_id).original_task_days_mask;
    muscimol_days_mask = passive_master_activity(animal_id).muscimol_days_mask;
    original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);
    first_day_idx = original_task_day_idx(1);
    
    % initialize avg_ROI trace to plot for this animal
    avg_this_ROI = nan(1,length(original_task_day_idx));
    
    for day_idx = original_task_day_idx
        
        % load activity
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        
        % skip if no activity
        if isempty(stim_act)
            continue
        end
        
        % deconvolve and baseline substract
        deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
        deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
        
        % get trial information
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        stim_wheel_move = passive_master_activity(animal_id).stim_wheel_move{day_idx};
        
        % activity for current stim and no move
        no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
        this_stim_act = deconvolved_stim_avg_act(:,:,no_move_trials&trialStimulusValue==possible_stimuli(stim_idx));
        
        % get avg ROI
        this_ROI = AP_svd_roi(U_master(:,:,1:num_comp),this_stim_act,[],[],roi.visual_left);
        this_ROI = permute(this_ROI,[2,3,1]);
        this_avg_ROI = mean(this_ROI(small_window_idx,:),2);
        avg_this_ROI(day_idx-first_day_idx+1) = mean(this_avg_ROI,1);
    end
    
    % to account for difference in training days
    try
        all_avg_visual_left(animal_id,:) = avg_this_ROI;
    catch
        all_avg_visual_left(animal_id,:) = [avg_this_ROI NaN];
    end      
end

original_mean_other_avg_visual_left = nanmean(all_avg_visual_left(2:end,:),1);
original_avg_visual_left = nanmean(all_avg_visual_left(1,:),1);
avg_visual_left = nanmean(all_avg_visual_left,1);

%% frontal left 

% stim left
stim_idx = 3;

% just initialize with nans for all animals and days of 107
% get day indexes
animal_id=1;
all_training_days = passive_master_activity(animal_id).day;
original_task_days_mask = passive_master_activity(animal_id).original_task_days_mask;
muscimol_days_mask = passive_master_activity(animal_id).muscimol_days_mask;
original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);
num_original_days = length(original_task_day_idx);

all_avg_frontal_left = nan(length(animals), length(original_task_day_idx));

% go through each animal
for animal_id=1:length(animals)
    
    animal = animals(animal_id);
    
    % get day indexes
    all_training_days = passive_master_activity(animal_id).day;
    original_task_days_mask = passive_master_activity(animal_id).original_task_days_mask;
    muscimol_days_mask = passive_master_activity(animal_id).muscimol_days_mask;
    original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);
    first_day_idx = original_task_day_idx(1);
    
    % initialize avg_ROI trace to plot for this animal
    avg_this_ROI = nan(1,length(original_task_day_idx));
    
    for day_idx = original_task_day_idx
        
        % load activity
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        
        % skip if no activity
        if isempty(stim_act)
            continue
        end
        
        % deconvolve and baseline substract
        deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
        deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
        
        % get trial information
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        stim_wheel_move = passive_master_activity(animal_id).stim_wheel_move{day_idx};
        
        % activity for current stim and no move
        no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
        this_stim_act = deconvolved_stim_avg_act(:,:,no_move_trials&trialStimulusValue==possible_stimuli(stim_idx));
        
        % get avg ROI
        this_ROI = AP_svd_roi(U_master(:,:,1:num_comp),this_stim_act,[],[],roi.frontal_left);
        this_ROI = permute(this_ROI,[2,3,1]);
        this_avg_ROI = mean(this_ROI(small_window_idx,:),2);
        avg_this_ROI(day_idx-first_day_idx+1) = mean(this_avg_ROI,1);
    end
    
    % to account for difference in training days
    try
        all_avg_frontal_left(animal_id,:) = avg_this_ROI;
    catch
        all_avg_frontal_left(animal_id,:) = [avg_this_ROI NaN];
    end      
end

original_mean_other_avg_frontal_left = nanmean(all_avg_frontal_left(2:end,:),1);
original_avg_frontal_left = nanmean(all_avg_frontal_left(1,:),1);
avg_frontal_left = nanmean(all_avg_frontal_left,1);

%% plot for average

%% - all mice
% visual
avg_visual_left = nanmean(all_avg_visual_left,1);
std_visual_left = std(all_avg_visual_left,1);
sem_visual_left = std_visual_left./sqrt(length(animals));

% colors
roi_colors = max(0, brewermap(1,'Oranges')-0.2);
% name of plot
figure('Name',['All avg fluorescence for visual left for stim ' num2str(possible_stimuli(stim_idx))]);
% plot
errorbar(1:num_original_days, avg_visual_left, sem_visual_left, '-o', 'LineWidth', 2);
hold on;
set(gca,'ColorOrder',roi_colors)
ylabel('ROI fluorescence');
xlabel('Training day');
ylim([-1*10^(-3) 4*10^(-3)]);
AM_prettyfig


% frontal
avg_frontal_left = nanmean(all_avg_frontal_left,1);
std_frontal_left = std(all_avg_frontal_left,1);
sem_frontal_left = std_frontal_left./sqrt(length(animals));

% colors
roi_colors = max(0, brewermap(1,'Reds')-0.4);
% name of plot
figure('Name',['All avg fluorescence for frontal left for stim ' num2str(possible_stimuli(stim_idx))]);
% plot
errorbar(1:num_original_days, avg_frontal_left, sem_frontal_left, '-o', 'LineWidth', 2);
hold on;
set(gca,'ColorOrder',roi_colors)
ylabel('ROI fluorescence');
xlabel('Training day');
ylim([-1*10^(-3) 4*10^(-3)]);
AM_prettyfig

%% - AP107
% colors
roi_colors = max(0, brewermap(1,'Oranges')-0.2);
% name of plot
figure('Name',['AP107 avg fluorescence for visual left for stim ' num2str(possible_stimuli(stim_idx))]);
% plot
plot(1:num_original_days, original_avg_visual_left, '-o', 'LineWidth', 2);
hold on;
set(gca,'ColorOrder',roi_colors)
ylabel('ROI fluorescence');
xlabel('Training day');
ylim([-1*10^(-3) 6*10^(-3)]);
AM_prettyfig

% AP107
% colors
roi_colors = max(0, brewermap(1,'Reds')-0.4);
% name of plot
figure('Name',['AP107 avg fluorescence for frontal left for stim ' num2str(possible_stimuli(stim_idx))]);
% plot
plot(1:num_original_days, original_avg_frontal_left, '-o', 'LineWidth', 2);
hold on;
set(gca,'ColorOrder',roi_colors)
ylabel('ROI fluorescence');
xlabel('Training day');
ylim([-1*10^(-3) 6*10^(-3)]);
AM_prettyfig

% %  ----- separated - by learner in reversal
% % colors
% colors_original_avg = cat(1, [0 0 0], [0.5 0.5 0.5]);
% % name of plot
% figure('Name',['ROI trace for visual left for stim ' num2str(possible_stimuli(stim_idx))]);
% ax = axes;
% % plot
% plot(1:num_original_days, original_avg_visual_left, '-o', 'LineWidth', 2);
% hold on;
% plot(1:num_original_days,original_mean_other_avg_visual_left, '-o', 'LineWidth', 2);
% hold on;
% ylabel('ROI fluorescence');
% xlabel('Training day');
% ylim([-2*10^(-3) 6*10^(-3)]);
% ax.ColorOrder = colors_original_avg;
% legend({'learner', 'non-learners'})
% AM_prettyfig

% 
% % ----- separated
% % colors
% colors_original_avg = cat(1, [0 0 0], [0.5 0.5 0.5]);
% % name of plot
% figure('Name',['ROI trace for frontal left for stim ' num2str(possible_stimuli(stim_idx))]);
% ax = axes;
% % plot
% plot(1:num_original_days, original_avg_frontal_left, '-o', 'LineWidth', 2);
% hold on;
% plot(1:num_original_days,original_mean_other_avg_frontal_left, '-o', 'LineWidth', 2);
% hold on;
% ylabel('ROI fluorescence');
% xlabel('Training day');
% ylim([-2*10^(-3) 6*10^(-3)]);
% ax.ColorOrder = colors_original_avg;
% legend({'learner', 'non-learners'})
% AM_prettyfig

%% STATS

[p,tbl,stats] = anova1(all_avg_visual_left(:,1:10));

[p,tbl,stats] = anova1(all_avg_frontal_left(:,1:10));

%% - Reversal 

%% visual right 
% stim left
stim_idx = 1;

% just initialize with nans for all animals and days of 115
% get day indexes
all_training_days = passive_master_activity(animal_id).day;
reversal_task_days_mask = passive_master_activity(animal_id).reversal_task_days_mask;
muscimol_days_mask = passive_master_activity(animal_id).muscimol_days_mask;
reversal_task_day_idx = find(reversal_task_days_mask-muscimol_days_mask==1);

all_avg_visual_right = nan(length(animals), length(reversal_task_day_idx));

% go through each animal
for animal_id=1:length(animals)
    
    animal = animals(animal_id);
    
    % get day indexes
    all_training_days = passive_master_activity(animal_id).day;
    reversal_task_days_mask = passive_master_activity(animal_id).reversal_task_days_mask;
    muscimol_days_mask = passive_master_activity(animal_id).muscimol_days_mask;
    reversal_task_day_idx = find(reversal_task_days_mask-muscimol_days_mask==1);
    first_day_idx = reversal_task_day_idx(1);
    
    % initialize avg_ROI trace to plot for this animal
    avg_this_ROI = nan(1,length(reversal_task_day_idx));
    
    for day_idx = reversal_task_day_idx
        
        % load activity
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        
        % skip if no activity
        if isempty(stim_act)
            continue
        end
        
        % deconvolve and baseline substract
        deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
        deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
        
        % get trial information
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        stim_wheel_move = passive_master_activity(animal_id).stim_wheel_move{day_idx};
        
        % activity for current stim and no move
        no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
        this_stim_act = deconvolved_stim_avg_act(:,:,no_move_trials&trialStimulusValue==possible_stimuli(stim_idx));
        
        % get avg ROI
        this_ROI = AP_svd_roi(U_master(:,:,1:num_comp),this_stim_act,[],[],roi.visual_right);
        this_ROI = permute(this_ROI,[2,3,1]);
        this_avg_ROI = mean(this_ROI(small_window_idx,:),2);
        avg_this_ROI(day_idx-first_day_idx+1) = mean(this_avg_ROI,1);
    end
    
    % to account for difference in training days
    try
        all_avg_visual_right(animal_id,:) = avg_this_ROI;
    catch
        all_avg_visual_right(animal_id,:) = [avg_this_ROI NaN NaN];
    end      
end

reversal_mean_other_avg_visual_right = nanmean(all_avg_visual_right(2:end,:),1);
reversal_avg_visual_right = nanmean(all_avg_visual_right(1,:),1);

%% frontal right 

% stim left
stim_idx = 1;

% just initialize with nans for all animals and days of 115
% get day indexes
all_training_days = passive_master_activity(animal_id).day;
reversal_task_days_mask = passive_master_activity(animal_id).reversal_task_days_mask;
muscimol_days_mask = passive_master_activity(animal_id).muscimol_days_mask;
reversal_task_day_idx = find(reversal_task_days_mask-muscimol_days_mask==1);

all_avg_frontal_right = nan(length(animals), length(reversal_task_day_idx));

% go through each animal
for animal_id=1:length(animals)
    
    animal = animals(animal_id);
    
    % get day indexes
    all_training_days = passive_master_activity(animal_id).day;
    reversal_task_days_mask = passive_master_activity(animal_id).reversal_task_days_mask;
    muscimol_days_mask = passive_master_activity(animal_id).muscimol_days_mask;
    reversal_task_day_idx = find(reversal_task_days_mask-muscimol_days_mask==1);
    first_day_idx = reversal_task_day_idx(1);
    
    % initialize avg_ROI trace to plot for this animal
    avg_this_ROI = nan(1,length(reversal_task_day_idx));
    
    for day_idx = reversal_task_day_idx
        
        % load activity
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        
        % skip if no activity
        if isempty(stim_act)
            continue
        end
        
        % deconvolve and baseline substract
        deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
        deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
        
        % get trial information
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        stim_wheel_move = passive_master_activity(animal_id).stim_wheel_move{day_idx};
        
        % activity for current stim and no move
        no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
        this_stim_act = deconvolved_stim_avg_act(:,:,no_move_trials&trialStimulusValue==possible_stimuli(stim_idx));
        
        % get avg ROI
        this_ROI = AP_svd_roi(U_master(:,:,1:num_comp),this_stim_act,[],[],roi.frontal_right);
        this_ROI = permute(this_ROI,[2,3,1]);
        this_avg_ROI = mean(this_ROI(small_window_idx,:),2);
        avg_this_ROI(day_idx-first_day_idx+1) = mean(this_avg_ROI,1);
    end
    
    % to account for difference in training days
    try
        all_avg_frontal_right(animal_id,:) = avg_this_ROI;
    catch
        all_avg_frontal_right(animal_id,:) = [avg_this_ROI NaN NaN];
    end      
end

reversal_mean_other_avg_frontal_right = nanmean(all_avg_frontal_right(2:end,:),1);
reversal_avg_frontal_right = nanmean(all_avg_frontal_right(1,:),1);

%% plot for average black and grey 

% visual
std_visual_right = std(all_avg_visual_right(2:end,:),1);
sem_visual_right = std_visual_right./sqrt(length(animals)-1);
% name of plot
figure('Name',['ROI trace for visual right for stim ' num2str(possible_stimuli(stim_idx))]);
% plot
plot(1:length(reversal_task_day_idx), reversal_avg_visual_right, '-or', 'LineWidth', 2);
hold on;
errorbar(1:length(reversal_task_day_idx),reversal_mean_other_avg_visual_right, sem_visual_right, '-ok', 'LineWidth', 2);
hold on;
ylabel('ROI fluorescence');
xlabel('Training day');
ylim([-2*10^(-3) 6*10^(-3)]);
legend({'learner', 'non-learners'})
AM_prettyfig


% frontal
std_frontal_right = std(all_avg_frontal_right(2:end,:),1);
sem_frontal_right = std_frontal_right./sqrt(length(animals)-1);
% name of plot
figure('Name',['ROI trace for frontal right for stim ' num2str(possible_stimuli(stim_idx))]);
% plot
plot(1:length(reversal_task_day_idx), reversal_avg_frontal_right, '-or', 'LineWidth', 2);
hold on;
errorbar(1:length(reversal_task_day_idx),reversal_mean_other_avg_frontal_right, sem_frontal_right, '-ok', 'LineWidth', 2);
hold on;
ylabel('ROI fluorescence');
xlabel('Training day');
ylim([-2*10^(-3) 6*10^(-3)]);
legend({'learner', 'non-learners'})
AM_prettyfig

% AP107
% colors
roi_colors = max(0, brewermap(1,'Blues')-0.5);
% name of plot
figure('Name',['AP107 avg fluorescence for frontal right for stim ' num2str(possible_stimuli(stim_idx))]);
% plot
plot(reversal_avg_frontal_right, '-o', 'LineWidth', 2);
hold on;
set(gca,'ColorOrder',roi_colors)
ylabel('ROI fluorescence');
xlabel('Training day');
ylim([-2*10^(-3) 6*10^(-3)]);
AM_prettyfig

% colors
roi_colors = max(0, brewermap(1,'Purples')-0.3);
% name of plot
figure('Name',['AP107 avg fluorescence for visual right for stim ' num2str(possible_stimuli(stim_idx))]);
% plot
plot(reversal_avg_visual_right, '-o', 'LineWidth', 2);
hold on;
set(gca,'ColorOrder',roi_colors)
ylabel('ROI fluorescence');
xlabel('Training day');
ylim([-2*10^(-3) 6*10^(-3)]);
AM_prettyfig

%% STATS

[p,tbl,stats] = anova1(all_avg_visual_right(2:end,1:10));

[p,tbl,stats] = anova1(all_avg_frontal_right(2:end,1:10));

% [p,tbl,stats] = anova1(all_avg_visual_right(1,1:10));


learner_frontal_right_mdl = fitlm(1:10, all_avg_frontal_right(1,1:10));
other_frontal_right_mdl = fitlm(1:10, reversal_mean_other_avg_frontal_right(1:10));

learner_visual_right_mdl = fitlm(1:10, all_avg_visual_right(1,1:10));
other_visual_right_mdl = fitlm(1:10, reversal_mean_other_avg_visual_right(1:10));


%% - Post-training and last day of training

% - visual left ROIs 
% get day indexes for 113
animal_id = 3;
all_training_days = passive_master_activity(animal_id).day;
post_training_passive_days_mask = passive_master_activity(animal_id).post_training_passive_days_mask;
post_training_passive_day_idx = find(post_training_passive_days_mask);

visual_left_avg_ROI = nan(length(post_training_passive_day_idx),length(timevec),length(animals));
for animal_id=3:length(animals)
    all_training_days = passive_master_activity(animal_id).day;
    post_training_passive_days_mask = passive_master_activity(animal_id).post_training_passive_days_mask;
    post_training_passive_day_idx = find(post_training_passive_days_mask);

    for day_idx=post_training_passive_day_idx
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
        deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
        
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        right_stim_act = deconvolved_stim_avg_act(:,:,trialStimulusValue==1);
        visual_left = AP_svd_roi(U_master(:,:,1:num_comp),right_stim_act,[],[],roi.visual_left);
        visual_left = permute(visual_left,[2,3,1]);
        visual_left_avg_ROI(day_idx-post_training_passive_day_idx(1)+1,:,animal_id) = mean(visual_left,2);
    end
end

% average across animals
this_visual_left_avg_ROI = visual_left_avg_ROI(:,:,3:end);
avg_visual_left_ROI = nanmean(this_visual_left_avg_ROI,3);
std_visual_left = std(this_visual_left_avg_ROI, 0, 3);
sem_visual_left =  std_visual_left./sqrt(length(animals)-2);

% get last day of training ROIs
last_day_training_idx = post_training_passive_day_idx(1)-1;

last_visual_left_ROI = nan(length(timevec),length(animals));
day_idx = last_day_training_idx;
for animal_id=3:length(animals)
    stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
    deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
    deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
    
    trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
    right_stim_act = deconvolved_stim_avg_act(:,:,trialStimulusValue==1);
    visual_left = AP_svd_roi(U_master(:,:,1:num_comp),right_stim_act,[],[],roi.visual_left);
    visual_left = permute(visual_left,[2,3,1]);
    last_visual_left_ROI(:,animal_id) = mean(visual_left,2);
end

% average across animals
this_visual_left_avg_ROI = last_visual_left_ROI(:,3:end);
avg_last_visual_left_ROI = nanmean(this_visual_left_avg_ROI,2);

% get all values
all_post_visual_left_ROI = [avg_last_visual_left_ROI(small_window_idx)'; avg_visual_left_ROI(:,small_window_idx)];
avg_all_post_visual_left = nanmean(all_post_visual_left_ROI, 2);
std_post_visual_left = std(all_post_visual_left_ROI, 0, 2);
sem_post_visual_left =  std_post_visual_left./sqrt(length(animals)-2);


% - frontal left ROIs 
% get day indexes for 113
animal_id = 3;
all_training_days = passive_master_activity(animal_id).day;
post_training_passive_days_mask = passive_master_activity(animal_id).post_training_passive_days_mask;
post_training_passive_day_idx = find(post_training_passive_days_mask);

frontal_left_avg_ROI = nan(length(post_training_passive_day_idx),length(timevec),length(animals));
for animal_id=3:length(animals)
    all_training_days = passive_master_activity(animal_id).day;
    post_training_passive_days_mask = passive_master_activity(animal_id).post_training_passive_days_mask;
    post_training_passive_day_idx = find(post_training_passive_days_mask);

    for day_idx=post_training_passive_day_idx
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
        deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
        
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        right_stim_act = deconvolved_stim_avg_act(:,:,trialStimulusValue==1);
        frontal_left = AP_svd_roi(U_master(:,:,1:num_comp),right_stim_act,[],[],roi.frontal_left);
        frontal_left = permute(frontal_left,[2,3,1]);
        frontal_left_avg_ROI(day_idx-post_training_passive_day_idx(1)+1,:,animal_id) = mean(frontal_left,2);
    end
end

% average across animals
this_frontal_left_avg_ROI = frontal_left_avg_ROI(:,:,3:end);
avg_frontal_left_ROI = nanmean(this_frontal_left_avg_ROI,3);
std_frontal_left = std(this_frontal_left_avg_ROI, 0, 3);
sem_frontal_left =  std_frontal_left./sqrt(length(animals)-2);


% get last day of training ROIs

last_day_training_idx = post_training_passive_day_idx(1)-1;

last_frontal_left_ROI = nan(length(timevec),length(animals));
day_idx = last_day_training_idx;
for animal_id=3:length(animals)
    stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
    deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
    deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
    
    trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
    right_stim_act = deconvolved_stim_avg_act(:,:,trialStimulusValue==1);
    frontal_left = AP_svd_roi(U_master(:,:,1:num_comp),right_stim_act,[],[],roi.frontal_left);
    frontal_left = permute(frontal_left,[2,3,1]);
    last_frontal_left_ROI(:,animal_id) = mean(frontal_left,2);
end

% average across animals
this_frontal_left_avg_ROI = last_frontal_left_ROI(:,3:end);
avg_last_frontal_left_ROI = nanmean(this_frontal_left_avg_ROI,2);

% get all values
all_post_frontal_left_ROI = [avg_last_frontal_left_ROI(small_window_idx)'; avg_frontal_left_ROI(:,small_window_idx)];
avg_all_post_frontal_left = nanmean(all_post_frontal_left_ROI, 2);
std_post_frontal_left = std(all_post_frontal_left_ROI, 0, 2);
sem_post_frontal_left =  std_post_frontal_left./sqrt(length(animals)-2);

%% STATS

% frontal left

tmp_last_frontal = nanmean(this_frontal_left_avg_ROI(small_window_idx,:),1);
tmp_post_frontal = permute(nanmean(frontal_left_avg_ROI(:,small_window_idx,3:end),2), [1, 3, 2]); % mouse is column!
tmp_all_frontal = [tmp_last_frontal; tmp_post_frontal];

[p,tbl,stats] = anova1(tmp_all_frontal');

% visual left
tmp_last_visual = nanmean(this_visual_left_avg_ROI(small_window_idx,:),1);
tmp_post_visual = permute(nanmean(visual_left_avg_ROI(:,small_window_idx,3:end),2), [1, 3, 2]); % mouse is column!
tmp_all_visual = [tmp_last_visual; tmp_post_visual];

[p,tbl,stats] = anova1(tmp_all_visual');


%% -- plot

% visual
roi_colors = max(0, brewermap(1,'Oranges')-0.2);
figure('Name',['Avg fluorescence post-train vis left for stim ' num2str(possible_stimuli(stim_idx))]);
errorbar(1:4, avg_all_post_visual_left, sem_post_visual_left, '-o', 'LineWidth', 2);
hold on;
set(gca,'ColorOrder',roi_colors)
ylabel('ROI fluorescence');
xlabel('Training day');
ylim([-2*10^(-3) 4*10^(-3)]);
xlim([0.5 4.5])
xticklabels({'', 'last day', '', 'week 1', '', 'week 2', '', 'week 3', ''})
AM_prettyfig

% frontal
roi_colors = max(0, brewermap(1,'Reds')-0.4);
figure('Name',['Avg fluorescence post-train vis left for stim ' num2str(possible_stimuli(stim_idx))]);
errorbar(1:4, avg_all_post_frontal_left, sem_post_frontal_left, '-o', 'LineWidth', 2);
hold on;
set(gca,'ColorOrder',roi_colors)
ylabel('ROI fluorescence');
xlabel('Training day');
ylim([-2*10^(-3) 4*10^(-3)]);
xlim([0.5 4.5])
xticklabels({'', 'last day', '', 'week 1', '', 'week 2', '', 'week 3', ''})
AM_prettyfig

%% plot for average


%% - all mice
% visual
avg_visual_left = nanmean(all_avg_visual_left,1);
std_visual_left = std(all_avg_visual_left,1);
sem_visual_left = std_visual_left./sqrt(length(animals));

% colors
roi_colors = max(0, brewermap(1,'Oranges')-0.2);
% name of plot
figure('Name',['All avg fluorescence for visual left for stim ' num2str(possible_stimuli(stim_idx))]);
% plot
errorbar(1:num_original_days, avg_visual_left, sem_visual_left, '-o', 'LineWidth', 2);
hold on;
set(gca,'ColorOrder',roi_colors)
ylabel('ROI fluorescence');
xlabel('Training day');
ylim([-1*10^(-3) 4*10^(-3)]);
AM_prettyfig


% frontal
avg_frontal_left = nanmean(all_avg_frontal_left,1);
std_frontal_left = std(all_avg_frontal_left,1);
sem_frontal_left = std_frontal_left./sqrt(length(animals));

% colors
roi_colors = max(0, brewermap(1,'Reds')-0.4);
% name of plot
figure('Name',['All avg fluorescence for frontal left for stim ' num2str(possible_stimuli(stim_idx))]);
% plot
errorbar(1:num_original_days, avg_frontal_left, sem_frontal_left, '-o', 'LineWidth', 2);
hold on;
set(gca,'ColorOrder',roi_colors)
ylabel('ROI fluorescence');
xlabel('Training day');
ylim([-1*10^(-3) 4*10^(-3)]);
AM_prettyfig



%% RT and Passive fluorescence ----------------------------------------------

% load activity from passive
load('passive_master_activity.mat');
animals = {passive_master_activity.animal};

timestep = passive_master_activity.timestep;
start_time = passive_master_activity.start_time;
end_time = passive_master_activity.end_time;
stim_frame = passive_master_activity.stim_frame;
timevec = passive_master_activity.timevec;
num_comp = passive_master_activity.num_comp;

% load ROI masks
load('ROIs.mat');
all_ROIs = fieldnames(roi);

% load master U
master_u_fn = 'D:\Andrada\Master project\widefield_alignment\U_master';
load(master_u_fn);

% choose averaging window
small_window = [0.05 0.2];
small_window_idx = timevec>=small_window(1)&timevec<=small_window(2);

% possible stims
possible_stimuli = [-1 0 1];
%% - original task

% stim right
stim_idx = 3;

% just initialize with nans for all animals and days of 107
% get day indexes
animal_id = 1;
all_training_days = passive_master_activity(animal_id).day;
original_task_days_mask = passive_master_activity(animal_id).original_task_days_mask;
muscimol_days_mask = passive_master_activity(animal_id).muscimol_days_mask;
original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);

all_avg_frontal_left = nan(length(animals), length(original_task_day_idx));

% go through each animal
for animal_id=1:length(animals)
    
    animal = animals(animal_id);
    
    % get day indexes
    all_training_days = passive_master_activity(animal_id).day;
    original_task_days_mask = passive_master_activity(animal_id).original_task_days_mask;
    muscimol_days_mask = passive_master_activity(animal_id).muscimol_days_mask;
    original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);
    first_day_idx = original_task_day_idx(1);
    
    % initialize avg_ROI trace to plot for this animal
    avg_this_ROI = nan(1,length(original_task_day_idx));
    
    for day_idx = original_task_day_idx
        
        % load activity
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        
        % skip if no activity
        if isempty(stim_act)
            continue
        end
        
        % deconvolve and baseline substract
        deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
        deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
        
        % get trial information
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        stim_wheel_move = passive_master_activity(animal_id).stim_wheel_move{day_idx};
        
        % activity for current stim and no move
        no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
        this_stim_act = deconvolved_stim_avg_act(:,:,no_move_trials&trialStimulusValue==possible_stimuli(stim_idx));
        
        % get avg ROI
        this_ROI = AP_svd_roi(U_master(:,:,1:num_comp),this_stim_act,[],[],roi.frontal_left);
        this_ROI = permute(this_ROI,[2,3,1]);
        this_avg_ROI = mean(this_ROI(small_window_idx,:),2);
        avg_this_ROI(day_idx-first_day_idx+1) = mean(this_avg_ROI,1);
    end
    
    % to account for difference in training days
    try
        all_avg_frontal_left(animal_id,:) = avg_this_ROI;
    catch 
        all_avg_frontal_left(animal_id,:) = [avg_this_ROI NaN];
    end      
end

original_mean_other_avg_frontal_left = nanmean(all_avg_frontal_left(2:end,:),1);
original_avg_frontal_left = nanmean(all_avg_frontal_left(1,:),1);

% visual 
% stim right
stim_idx = 3;

% just initialize with nans for all animals and days of 107
% get day indexes
animal_id = 1;
all_training_days = passive_master_activity(animal_id).day;
original_task_days_mask = passive_master_activity(animal_id).original_task_days_mask;
muscimol_days_mask = passive_master_activity(animal_id).muscimol_days_mask;
original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);

all_avg_visual_left = nan(length(animals), length(original_task_day_idx));

% go through each animal
for animal_id=1:length(animals)
    
    animal = animals(animal_id);
    
    % get day indexes
    all_training_days = passive_master_activity(animal_id).day;
    original_task_days_mask = passive_master_activity(animal_id).original_task_days_mask;
    muscimol_days_mask = passive_master_activity(animal_id).muscimol_days_mask;
    original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);
    first_day_idx = original_task_day_idx(1);
    
    % initialize avg_ROI trace to plot for this animal
    avg_this_ROI = nan(1,length(original_task_day_idx));
    
    for day_idx = original_task_day_idx
        
        % load activity
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        
        % skip if no activity
        if isempty(stim_act)
            continue
        end
        
        % deconvolve and baseline substract
        deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
        deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
        
        % get trial information
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        stim_wheel_move = passive_master_activity(animal_id).stim_wheel_move{day_idx};
        
        % activity for current stim and no move
        no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
        this_stim_act = deconvolved_stim_avg_act(:,:,no_move_trials&trialStimulusValue==possible_stimuli(stim_idx));
        
        % get avg ROI
        this_ROI = AP_svd_roi(U_master(:,:,1:num_comp),this_stim_act,[],[],roi.visual_left);
        this_ROI = permute(this_ROI,[2,3,1]);
        this_avg_ROI = mean(this_ROI(small_window_idx,:),2);
        avg_this_ROI(day_idx-first_day_idx+1) = mean(this_avg_ROI,1);
    end
    
    % to account for difference in training days
    try
        all_avg_visual_left(animal_id,:) = avg_this_ROI;
    catch 
        all_avg_visual_left(animal_id,:) = [avg_this_ROI NaN];
    end      
end

original_mean_other_avg_visual_left = nanmean(all_avg_visual_left(2:end,:),1);
original_avg_visual_left = nanmean(all_avg_visual_left(1,:),1);

%% - reversal task

% stim left
stim_idx = 1;

% just initialize with nans for all animals and days of 115
% get day indexes
animal_id = 5;
all_training_days = passive_master_activity(animal_id).day;
reversal_task_days_mask = passive_master_activity(animal_id).reversal_task_days_mask;
muscimol_days_mask = passive_master_activity(animal_id).muscimol_days_mask;
reversal_task_day_idx = find(reversal_task_days_mask-muscimol_days_mask==1);

all_avg_frontal_right = nan(length(animals), length(reversal_task_day_idx));

% go through each animal
for animal_id=1:length(animals)
    
    animal = animals(animal_id);
    
    % get day indexes
    all_training_days = passive_master_activity(animal_id).day;
    reversal_task_days_mask = passive_master_activity(animal_id).reversal_task_days_mask;
    muscimol_days_mask = passive_master_activity(animal_id).muscimol_days_mask;
    reversal_task_day_idx = find(reversal_task_days_mask-muscimol_days_mask==1);
    first_day_idx = reversal_task_day_idx(1);
    
    % initialize avg_ROI trace to plot for this animal
    avg_this_ROI = nan(1,length(reversal_task_day_idx));
    
    for day_idx = reversal_task_day_idx
        
        % load activity
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        
        % skip if no activity
        if isempty(stim_act)
            continue
        end
        
        % deconvolve and baseline substract
        deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
        deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
        
        % get trial information
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        stim_wheel_move = passive_master_activity(animal_id).stim_wheel_move{day_idx};
        
        % activity for current stim and no move
        no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
        this_stim_act = deconvolved_stim_avg_act(:,:,no_move_trials&trialStimulusValue==possible_stimuli(stim_idx));
        
        % get avg ROI
        this_ROI = AP_svd_roi(U_master(:,:,1:num_comp),this_stim_act,[],[],roi.frontal_right);
        this_ROI = permute(this_ROI,[2,3,1]);
        this_avg_ROI = mean(this_ROI(small_window_idx,:),2);
        avg_this_ROI(day_idx-first_day_idx+1) = mean(this_avg_ROI,1);
    end
    
    % to account for difference in training days
    try
        all_avg_frontal_right(animal_id,:) = avg_this_ROI;
    catch
        all_avg_frontal_right(animal_id,:) = [avg_this_ROI NaN NaN];
    end      
end

reversal_mean_other_avg_frontal_right = nanmean(all_avg_frontal_right(2:end,:),1);
reversal_avg_frontal_right = nanmean(all_avg_frontal_right(1,:),1);

% visual

% stim right
stim_idx = 1;

% just initialize with nans for all animals and days of 115
% get day indexes
animal_id = 5;
all_training_days = passive_master_activity(animal_id).day;
reversal_task_days_mask = passive_master_activity(animal_id).reversal_task_days_mask;
muscimol_days_mask = passive_master_activity(animal_id).muscimol_days_mask;
reversal_task_day_idx = find(reversal_task_days_mask-muscimol_days_mask==1);

all_avg_visual_right = nan(length(animals), length(reversal_task_day_idx));

% go through each animal
for animal_id=1:length(animals)
    
    animal = animals(animal_id);
    
    % get day indexes
    all_training_days = passive_master_activity(animal_id).day;
    reversal_task_days_mask = passive_master_activity(animal_id).reversal_task_days_mask;
    muscimol_days_mask = passive_master_activity(animal_id).muscimol_days_mask;
    reversal_task_day_idx = find(reversal_task_days_mask-muscimol_days_mask==1);
    first_day_idx = reversal_task_day_idx(1);
    
    % initialize avg_ROI trace to plot for this animal
    avg_this_ROI = nan(1,length(reversal_task_day_idx));
    
    for day_idx = reversal_task_day_idx
        
        % load activity
        stim_act = passive_master_activity(animal_id).stim_activity{day_idx};
        
        % skip if no activity
        if isempty(stim_act)
            continue
        end
        
        % deconvolve and baseline substract
        deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
        deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
        
        % get trial information
        trialStimulusValue = passive_master_activity(animal_id).trial_id{day_idx};
        stim_wheel_move = passive_master_activity(animal_id).stim_wheel_move{day_idx};
        
        % activity for current stim and no move
        no_move_trials = sum(stim_wheel_move(stim_frame:end,:),1)==0;
        this_stim_act = deconvolved_stim_avg_act(:,:,no_move_trials&trialStimulusValue==possible_stimuli(stim_idx));
        
        % get avg ROI
        this_ROI = AP_svd_roi(U_master(:,:,1:num_comp),this_stim_act,[],[],roi.visual_right);
        this_ROI = permute(this_ROI,[2,3,1]);
        this_avg_ROI = mean(this_ROI(small_window_idx,:),2);
        avg_this_ROI(day_idx-first_day_idx+1) = mean(this_avg_ROI,1);
    end
    
    % to account for difference in training days
    try
        all_avg_visual_right(animal_id,:) = avg_this_ROI;
    catch
        all_avg_visual_right(animal_id,:) = [avg_this_ROI NaN NaN];
    end      
end

reversal_mean_other_avg_visual_right = nanmean(all_avg_visual_right(2:end,:),1);
reversal_avg_visual_right = nanmean(all_avg_visual_right(1,:),1);


%% - RT
load('all_mice_behaviour.mat');

% Original task 
animal_id = 1;
all_training_days = behaviour(1).day;
original_task_days_mask = behaviour(1).original_task_days_mask;
muscimol_days =  behaviour(1).muscimol_days;
muscimol_days_mask = ismember(all_training_days, muscimol_days);
original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);

original_all_RT = nan(length(animals), length(original_task_day_idx));
for animal_id=1:length(animals)
    all_training_days = behaviour(animal_id).day;
    original_task_days_mask = behaviour(animal_id).original_task_days_mask;
    muscimol_days =  behaviour(animal_id).muscimol_days;
    muscimol_days_mask = ismember(all_training_days, muscimol_days);
    original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);
    
    % get reaction times percentage that falls between 100-200 ms
    reaction_times_percentage = nan(1,length(original_task_day_idx));
    first_day_idx = original_task_day_idx(1);
    for day_idx = original_task_day_idx
        reaction_times = behaviour(animal_id).reaction_times{day_idx};
        reaction_times_percentage(day_idx-first_day_idx+1) = sum(discretize(reaction_times,[0.1 0.2])==1)/length(reaction_times);
    end
    try
        original_all_RT(animal_id, :) = reaction_times_percentage;
    catch
        original_all_RT(animal_id, :) = [reaction_times_percentage NaN];
    end
end

% ger average RT traces
original_other_mean_RT = nanmean(original_all_RT(2:end,:), 1);
original_RT = nanmean(original_all_RT(1,:), 1);

% Reversal 
animal_id = 5;

% days for 115
all_training_days = behaviour(animal_id).day;
reversal_task_days_mask = behaviour(animal_id).reversal_task_days_mask;
muscimol_days =  behaviour(animal_id).muscimol_days;
muscimol_days_mask = ismember(all_training_days, muscimol_days);
reversal_task_day_idx = find(reversal_task_days_mask-muscimol_days_mask==1);

reversal_all_RT = nan(length(animals), length(reversal_task_day_idx));
for animal_id=1:length(animals)
    all_training_days = behaviour(animal_id).day;
    reversal_task_days_mask = behaviour(animal_id).reversal_task_days_mask;
    muscimol_days =  behaviour(animal_id).muscimol_days;
    muscimol_days_mask = ismember(all_training_days, muscimol_days);
    reversal_task_day_idx = find(reversal_task_days_mask-muscimol_days_mask==1);
    
    % get reaction times percentage that falls between 100-200 ms
    reaction_times_percentage = nan(1,length(reversal_task_day_idx));
    first_day_idx = reversal_task_day_idx(1);
    for day_idx = reversal_task_day_idx
        reaction_times = behaviour(animal_id).reaction_times{day_idx};
        reaction_times_percentage(day_idx-first_day_idx+1) = sum(discretize(reaction_times,[0.1 0.2])==1)/length(reaction_times);
    end
    try
        reversal_all_RT(animal_id, :) = reaction_times_percentage;
    catch
        reversal_all_RT(animal_id, :) = [reaction_times_percentage NaN NaN];
    end
end

% ger average RT traces
reversal_other_mean_RT = nanmean(reversal_all_RT(2:end,:), 1);
reversal_RT = nanmean(reversal_all_RT(1,:), 1);

%% - plot

% original
roi_colors = max(0, brewermap(2,'Reds')-0.3);
figure;
title('Original RT and fluorescence')
hold on;
% RT
yyaxis left
set(gca,'ColorOrder',roi_colors)
plot(original_RT);
hold on;
plot(original_other_mean_RT)
hold on;
hline = findobj(gcf, 'type', 'line');
set(hline(1),'LineStyle','--')
ylim([0 0.8])
ylabel('Percentage of reaction times within 100-200 ms')
% fluorescence
yyaxis right
set(gca,'ColorOrder',roi_colors)
plot(original_avg_frontal_left)
hold on;
plot(original_mean_other_avg_frontal_left)
hold on;
hline = findobj(gcf, 'type', 'line');
set(hline(2),'LineStyle','--')
ylim([-1*10^-3 3*10^-3])
ylabel('Average fluorescence across training')
% make axes black
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
% legend
legend({'learner RT', 'learner fluorescence', 'non-learners RT', 'non-learners fluorescence'})
AM_prettyfig


% reversal
roi_colors = max(0, brewermap(2,'Blues')-0.3);
figure;
title('Reversal RT and fluorescence')
hold on;
% RT
yyaxis left
set(gca,'ColorOrder',roi_colors)
plot(reversal_RT);
hold on;
plot(reversal_other_mean_RT)
hold on;
hline = findobj(gcf, 'type', 'line');
set(hline(1),'LineStyle','--')
ylim([0 0.8])
ylabel('Percentage of reaction times within 100-200 ms')
% fluorescence
yyaxis right
set(gca,'ColorOrder',roi_colors)
plot(reversal_avg_frontal_right)
hold on;
plot(reversal_mean_other_avg_frontal_right)
hold on;
hline = findobj(gcf, 'type', 'line');
set(hline(2),'LineStyle','--')
ylim([-1*10^-3 3*10^-3])
ylabel('Average fluorescence across training')
% make axes black
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
% legend
legend({'learner RT', 'learner fluorescence', 'non-learners RT', 'non-learners fluorescence'})
AM_prettyfig


% original - make different plot
colors_learner = brewermap(length(original_other_mean_RT), 'Reds');
figure;
scatter(original_mean_other_avg_frontal_left, original_other_mean_RT, 100, colors_learner, 'filled')
hold on;
scatter(original_avg_frontal_left, original_RT, 100, colors_learner, 'filled')
xlabel('Fluorescence')
ylabel('Percentage of RT in 100-200 ms')
lsline
AM_prettyfig


% reversal - make different plot
colors_learner = brewermap(length(reversal_other_mean_RT), 'Reds');
colors_non_learners = brewermap(length(reversal_other_mean_RT), 'Greys');
figure;
scatter(reversal_mean_other_avg_frontal_right, reversal_other_mean_RT, 100, colors_non_learners, 'filled')
lsline
hold on;
scatter(reversal_avg_frontal_right, reversal_RT, 100, colors_learner, 'filled')
lsline
xlabel('Fluorescence')
ylabel('Percentage of RT in 100-200 ms')
lsline
AM_prettyfig

% NEW ---------------------------------------------------------------------
% FRONTAL
% original - all mice together
original_all_avg_frontal_left = mean([original_avg_frontal_left; original_mean_other_avg_frontal_left]);
original_all_RT = mean([original_RT; original_other_mean_RT]);

colors_learner = brewermap(length(original_other_mean_RT), 'Reds');
figure;
scatter(original_all_avg_frontal_left, original_all_RT, 100, colors_learner, 'filled')
h1 = lsline;
h1.Color = 'r';
xlabel('Fluorescence')
ylabel('Percentage of RT in 100-200 ms')
ylim([0 0.7])
xlim([-0.001 0.035])
AM_prettyfig

% reversal - make separate plots
colors_learner = brewermap(length(reversal_other_mean_RT), 'Reds');
colors_non_learners = brewermap(length(reversal_other_mean_RT), 'Greys');

figure;
scatter(reversal_mean_other_avg_frontal_right, reversal_other_mean_RT, 100, colors_non_learners, 'filled')
h1 = lsline;
h1.Color = 'k';
xlabel('Fluorescence')
ylabel('Percentage of RT in 100-200 ms')
ylim([0 0.7])
xlim([-0.001 0.002])
AM_prettyfig

figure;
scatter(reversal_avg_frontal_right, reversal_RT, 100, colors_learner, 'filled')
h1 = lsline;
h1.Color = 'r';
xlabel('Fluorescence')
ylabel('Percentage of RT in 100-200 ms')
ylim([0 0.7])
xlim([-0.001 0.002])
AM_prettyfig

% VISUAL
% original - all mice together
original_all_avg_visual_left = mean([original_avg_visual_left; original_mean_other_avg_visual_left]);
original_all_RT = mean([original_RT; original_other_mean_RT]);

colors_learner = brewermap(length(original_other_mean_RT), 'Reds');
figure;
scatter(original_all_avg_visual_left, original_all_RT, 100, colors_learner, 'filled')
h1 = lsline;
h1.Color = 'r';
xlabel('Fluorescence')
ylabel('Percentage of RT in 100-200 ms')
ylim([0 0.7])
% xlim([-0.001 0.0035])
AM_prettyfig

% reversal - make separate plots
colors_learner = brewermap(length(reversal_other_mean_RT), 'Reds');
colors_non_learners = brewermap(length(reversal_other_mean_RT), 'Greys');

figure;
scatter(reversal_mean_other_avg_visual_right, reversal_other_mean_RT, 100, colors_non_learners, 'filled')
h1 = lsline;
h1.Color = 'k';
xlabel('Fluorescence')
ylabel('Percentage of RT in 100-200 ms')
ylim([0 0.7])
xlim([0 0.0035])
AM_prettyfig

figure;
scatter(reversal_avg_visual_right, reversal_RT, 100, colors_learner, 'filled')
h1 = lsline;
h1.Color = 'r';
xlabel('Fluorescence')
ylabel('Percentage of RT in 100-200 ms')
ylim([0 0.7])
xlim([0 0.0035])
AM_prettyfig

%% STATS

original_all_avg_visual_left = mean([original_avg_visual_left; original_mean_other_avg_visual_left]);
original_all_RT = mean([original_RT; original_other_mean_RT]);

original_all_avg_frontal_left = mean([original_avg_frontal_left; original_mean_other_avg_frontal_left]);
original_all_RT = mean([original_RT; original_other_mean_RT]);

% visual
mdl_reversal_visual_RT = fitlm(reversal_avg_visual_right, reversal_RT);

mdl_reversal_visual_other_RT = fitlm(reversal_mean_other_avg_visual_right, reversal_other_mean_RT);

mdl_original_visual_RT = fitlm(original_all_avg_visual_left, original_all_RT);

% frontal
mdl_reversal_frontal_RT = fitlm(reversal_avg_frontal_right, reversal_RT);

mdl_reversal_frontal_other_RT = fitlm(reversal_mean_other_avg_frontal_right, reversal_other_mean_RT);

mdl_original_frontal_RT = fitlm(original_all_avg_frontal_left, original_all_RT);



%% Average image with atlas
animal = 'AP107';
day = '2021-12-17';
experiment = 2; 

% load experiment
load_parts.imaging = true;
load_parts.cam = true;
verbose = true;
AP_load_experiment;

aligned_avg_im = AP_align_widefield(avg_im,animal,day);

figure; 
imagesc(aligned_avg_im); 
colormap(gray); 
axis image; 
axis off; 
hold on;    
% AP_reference_outline('ccf_aligned','#808080');

AP_reference_outline('ccf_aligned_lefthemi','#808080');