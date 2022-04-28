%% TASK FLUORESCENCE --------------------------------------------------------------

%% To see what ROI window to use 
% AP107
% load activity from task
load('task_master_activity.mat');
animals = {task_master_activity.animal};

timestep = task_master_activity.timestep;
start_time = task_master_activity.start_time;
end_time = task_master_activity.end_time;
stim_frame = task_master_activity.stim_frame;
timevec = task_master_activity.timevec;
num_comp = task_master_activity.num_comp;

% AP107
animal_id = 1;

% load ROI masks
load('ROIs.mat');

% load master U
master_u_fn = 'D:\Andrada\Master project\widefield_alignment\U_master';
load(master_u_fn);

% frontal left ROIs for each day

% get day indexes for original
all_training_days = task_master_activity(animal_id).day;
original_task_days_mask = task_master_activity(animal_id).original_task_days_mask;
muscimol_days_mask = task_master_activity(animal_id).muscimol_days_mask;
original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);

frontal_left_avg_ROI = nan(length(original_task_day_idx),length(timevec));
for day_idx=original_task_day_idx
    stim_act = task_master_activity(animal_id).stim_activity{day_idx};
    deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
    deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
    
    trialStimulusValue = task_master_activity(animal_id).trial_id{day_idx};
    % only completed trials
    trialStimulusValue = trialStimulusValue(1:size(stim_act,3));
    
    right_stim_act = deconvolved_stim_avg_act(:,:,trialStimulusValue==1);
    frontal_left = AP_svd_roi(U_master(:,:,1:num_comp),right_stim_act,[],[],roi.frontal_left);
    frontal_left = permute(frontal_left,[2,3,1]);
    frontal_left_avg_ROI(day_idx-original_task_day_idx(1)+1,:) = mean(frontal_left,2);   
end

roi_colors = brewermap(length(original_task_day_idx),'Reds');
figure; title('Task original each day left ROI');
hold on;
set(gca,'ColorOrder',roi_colors)
plot(timevec,frontal_left_avg_ROI);

% frontal right ROIs

% get day indexes
all_training_days = task_master_activity(animal_id).day;
reversal_task_days_mask = task_master_activity(animal_id).reversal_task_days_mask;
reversal_task_day_idx = find(reversal_task_days_mask);

frontal_right_avg_ROI = nan(length(reversal_task_day_idx),length(timevec));
for day_idx=reversal_task_day_idx
    stim_act = task_master_activity(animal_id).stim_activity{day_idx};
    deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
    deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
    
    trialStimulusValue = task_master_activity(animal_id).trial_id{day_idx};
    % only completed trials
    trialStimulusValue = trialStimulusValue(1:size(stim_act,3));
    
    right_stim_act = deconvolved_stim_avg_act(:,:,trialStimulusValue==-1);
    frontal_right = AP_svd_roi(U_master(:,:,1:num_comp),right_stim_act,[],[],roi.frontal_right);
    frontal_right = permute(frontal_right,[2,3,1]);
    frontal_right_avg_ROI(day_idx-reversal_task_day_idx(1)+1,:) = mean(frontal_right,2);   
end

roi_colors = brewermap(length(reversal_task_day_idx),'Blues');
figure; title('Task reversal each day right ROI');
hold on;
set(gca,'ColorOrder',roi_colors)
plot(timevec,frontal_right_avg_ROI);


%% Avg fluorescence
% averaging window
small_window = [0.05 0.3];
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
title('AP107 mPFC avg task fluorescence')
xlabel('Training day')
ylabel('Avg fluorescence 5-300 ms after stim onset')
legend({'Original', 'Reversal'})

%% - Avg fluorescence for both side ROIs

% AP107
animal_id = 1;

% window
small_window = [0.05 0.2];
small_window_idx = timevec>=small_window(1)&timevec<=small_window(2);

%% -- original

% get day indexes for original
all_training_days = task_master_activity(animal_id).day;
original_task_days_mask = task_master_activity(animal_id).original_task_days_mask;
muscimol_days_mask = task_master_activity(animal_id).muscimol_days_mask;
original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);

% initialize avg ROIs and avg fluorescence
original_frontal_left_avg_ROI = nan(length(original_task_day_idx),length(timevec));
original_frontal_right_avg_ROI = nan(length(original_task_day_idx),length(timevec));

avg_original_frontal_right = nan(1,length(original_task_day_idx));
avg_original_frontal_left = nan(1,length(original_task_day_idx));

for day_idx=original_task_day_idx
    stim_act = task_master_activity(animal_id).stim_activity{day_idx};
    deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
    deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
    
    trialStimulusValue = task_master_activity(animal_id).trial_id{day_idx};
    % only completed trials
    trialStimulusValue = trialStimulusValue(1:size(stim_act,3));
    
    right_stim_act = deconvolved_stim_avg_act(:,:,trialStimulusValue==1);
    frontal_left = AP_svd_roi(U_master(:,:,1:num_comp),right_stim_act,[],[],roi.frontal_left);
    frontal_left = permute(frontal_left,[2,3,1]);
    
    right_stim_act = deconvolved_stim_avg_act(:,:,trialStimulusValue==1);
    frontal_right = AP_svd_roi(U_master(:,:,1:num_comp),right_stim_act,[],[],roi.frontal_right);
    frontal_right = permute(frontal_right,[2,3,1]);
    
    % save this day for left/right ROIs and avg fluorescence 
    this_idx = day_idx-original_task_day_idx(1)+1;
    
    original_frontal_left_avg_ROI(this_idx,:) = mean(frontal_left,2);
    avg_original_frontal_left(this_idx) = mean(original_frontal_left_avg_ROI(this_idx,small_window_idx),2);

    original_frontal_right_avg_ROI(this_idx,:) = mean(frontal_right,2);
    avg_original_frontal_right(this_idx) = mean(original_frontal_right_avg_ROI(this_idx,small_window_idx),2);
end

% plot
figure;
plot(1:length(original_task_day_idx), avg_original_frontal_left, 'r-o')
hold on;
plot(1:length(original_task_day_idx), avg_original_frontal_right, 'b-o')
title('AP107 mPFC avg original task fluorescence')
xlabel('Training day')
ylabel('Avg fluorescence 5-200 ms after stim onset')
legend({'Left', 'Right'})
ylim([-1*10^-3 6*10^-3])


%% -- reversal

% get day indexes
all_training_days = task_master_activity(animal_id).day;
reversal_task_days_mask = task_master_activity(animal_id).reversal_task_days_mask;
reversal_task_day_idx = find(reversal_task_days_mask);

% initialize avg ROIs and avg fluorescence
reversal_frontal_left_avg_ROI = nan(length(reversal_task_day_idx),length(timevec));
reversal_frontal_right_avg_ROI = nan(length(reversal_task_day_idx),length(timevec));

avg_reversal_frontal_right = nan(1,length(reversal_task_day_idx));
avg_reversal_frontal_left = nan(1,length(reversal_task_day_idx));

for day_idx=reversal_task_day_idx
    stim_act = task_master_activity(animal_id).stim_activity{day_idx};
    deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
    deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
    
    trialStimulusValue = task_master_activity(animal_id).trial_id{day_idx};
    % only completed trials
    trialStimulusValue = trialStimulusValue(1:size(stim_act,3));
    
   	left_stim_act = deconvolved_stim_avg_act(:,:,trialStimulusValue==-1);
    frontal_left = AP_svd_roi(U_master(:,:,1:num_comp),left_stim_act,[],[],roi.frontal_left);
    frontal_left = permute(frontal_left,[2,3,1]);
    
    left_stim_act = deconvolved_stim_avg_act(:,:,trialStimulusValue==-1);
    frontal_right = AP_svd_roi(U_master(:,:,1:num_comp),left_stim_act,[],[],roi.frontal_right);
    frontal_right = permute(frontal_right,[2,3,1]);
    
    % save this day for left/right ROIs and avg fluorescence 
    this_idx = day_idx-reversal_task_day_idx(1)+1;
    
    reversal_frontal_left_avg_ROI(this_idx,:) = mean(frontal_left,2);
    avg_reversal_frontal_left(this_idx) = mean(reversal_frontal_left_avg_ROI(this_idx,small_window_idx),2);

    reversal_frontal_right_avg_ROI(this_idx,:) = mean(frontal_right,2);
    avg_reversal_frontal_right(this_idx) = mean(reversal_frontal_right_avg_ROI(this_idx,small_window_idx),2);
end

% plot
figure;
plot(1:length(reversal_task_day_idx), avg_reversal_frontal_left, 'r-o')
hold on;
plot(1:length(reversal_task_day_idx), avg_reversal_frontal_right, 'b-o')
title('AP107 mPFC avg reversal task fluorescence')
xlabel('Training day')
ylabel('Avg fluorescence 5-200 ms after stim onset')
legend({'Left', 'Right'})
ylim([-1*10^-3 6*10^-3])


%% Avg picture with scroll through days

% load activity from task
load('task_master_activity.mat');
animals = {task_master_activity.animal};

timestep = task_master_activity.timestep;
start_time = task_master_activity.start_time;
end_time = task_master_activity.end_time;
stim_frame = task_master_activity.stim_frame;
timevec = task_master_activity.timevec;
num_comp = task_master_activity.num_comp;

% choose averaging window
small_window = [0.05 0.3];
small_window_idx = timevec>=small_window(1)&timevec<=small_window(2);

% load master U
master_u_fn = 'D:\Andrada\Master project\widefield_alignment\U_master';
load(master_u_fn);

% possible stims
possible_stimuli = [-1 0 1];

for animal_id=1:length(animals)
    training_days = task_master_activity(animal_id).day;
    muscimol_days = training_days(task_master_activity(animal_id).muscimol_days_mask);
        
    % initialize avg_act over window
    all_window_avg_act = nan(num_comp, length(training_days));
    
    for day_idx=1:length(training_days)
        day = training_days(day_idx);
        
        % skip if it's a muscimol day
        is_muscimol = find(contains(muscimol_days,day));
        if is_muscimol
            continue
        end
        
        % get stim act for this animal and day
        stim_act = task_master_activity(animal_id).stim_activity{day_idx};
        
        % load trial information and wheel move
        trialStimulusValue = task_master_activity(animal_id).trial_id{day_idx};
        % only completed trials
        trialStimulusValue = trialStimulusValue(1:size(stim_act,3));
        
        stim_wheel_move = task_master_activity(animal_id).stim_wheel_move{day_idx};
        
        % find average activity for each stimulus
        stim_avg_act = nanmean(stim_act,3);
        stim_avg_act = stim_avg_act - stim_avg_act(:,stim_frame,:);
        
        % deconvolve
        deconvolved_stim_avg_act = AP_deconv_wf(stim_avg_act, [], 1/timestep);
        
        % average across chosen window
        window_avg_act = mean(deconvolved_stim_avg_act(:,small_window_idx,:),2);
        all_window_avg_act(:,day_idx,:) = window_avg_act;
    end
    
    % get fluoresence
    all_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),all_window_avg_act);
    %task_master_activity(animal_id).all_avg_fluorescence = all_avg_fluorescence;
    
    % video
    AP_image_scroll(all_avg_fluorescence,1:length(task_master_activity(animal_id).day)); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)])
    axis image;
end

% save the avg_fluorescence pics in struct to load easily
%save('task_master_activity.mat', 'task_master_activity', '-v7.3')

%% Frontal ROIs avg post-learning per task

figsdir = 'D:\Andrada\Master project\Random figures\ROI plots per task\';

% load activity from task
load('task_master_activity.mat');
animals = {task_master_activity.animal};

timestep = task_master_activity.timestep;
start_time = task_master_activity.start_time;
end_time = task_master_activity.end_time;
stim_frame = task_master_activity.stim_frame;
timevec = task_master_activity.timevec;
num_comp = task_master_activity.num_comp;

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
original_left_frontal_ROIs = nan(length(animals),length(timevec));
original_right_frontal_ROIs = nan(length(animals),length(timevec));

for animal_id=1:length(animals)
    
    % get day indexes
    all_training_days = task_master_activity(animal_id).day;
    original_task_days_mask = task_master_activity(animal_id).original_task_days_mask;
    muscimol_days_mask = task_master_activity(animal_id).muscimol_days_mask;
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
        stim_act = task_master_activity(animal_id).stim_activity{day_idx};
        
        % deconvolve and baseline substract
        deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
        deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
        
        % get avg left ROI for left stim
        left_ROI = AP_svd_roi(U_master(:,:,1:num_comp),deconvolved_stim_avg_act,[],[],roi.frontal_left);
        left_ROI = permute(left_ROI,[2,3,1]);
        this_left_ROI(:,day_idx-these_day_idx(1)+1) = mean(left_ROI,2);
        
        % get avg right ROI for left stim
        right_ROI = AP_svd_roi(U_master(:,:,1:num_comp),deconvolved_stim_avg_act,[],[],roi.frontal_right);
        right_ROI = permute(right_ROI,[2,3,1]);
        this_right_ROI(:,day_idx-these_day_idx(1)+1) = mean(right_ROI,2);
    end
    
    original_left_frontal_ROIs(animal_id,:) = mean(this_left_ROI,2);
    original_right_frontal_ROIs(animal_id,:) = mean(this_right_ROI,2);
end

%% - reversal task

% define learned day
learned_day_idx = 7;

% initialize frontal ROIs
reversal_left_frontal_ROIs = nan(length(animals),length(timevec));
reversal_right_frontal_ROIs = nan(length(animals),length(timevec));

for animal_id=1:length(animals)
    
    % get day indexes
    all_training_days = task_master_activity(animal_id).day;
    reversal_task_days_mask = task_master_activity(animal_id).reversal_task_days_mask;
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
        stim_act = task_master_activity(animal_id).stim_activity{day_idx};
        
        % deconvolve and baseline substract
        deconvolved_stim_avg_act = AP_deconv_wf(stim_act, [], 1/timestep);
        deconvolved_stim_avg_act = deconvolved_stim_avg_act - deconvolved_stim_avg_act(:,stim_frame,:);
        
        % get avg left ROI for left stim
        left_ROI = AP_svd_roi(U_master(:,:,1:num_comp),deconvolved_stim_avg_act,[],[],roi.frontal_left);
        left_ROI = permute(left_ROI,[2,3,1]);
        this_left_ROI(:,day_idx-these_day_idx(1)+1) = mean(left_ROI,2);
        
        % get avg right ROI for left stim
        right_ROI = AP_svd_roi(U_master(:,:,1:num_comp),deconvolved_stim_avg_act,[],[],roi.frontal_right);
        right_ROI = permute(right_ROI,[2,3,1]);
        this_right_ROI(:,day_idx-these_day_idx(1)+1) = mean(right_ROI,2);
    end
    
    reversal_left_frontal_ROIs(animal_id,:) = mean(this_left_ROI,2);
    reversal_right_frontal_ROIs(animal_id,:) = mean(this_right_ROI,2);
end

%% - plot for learner

figure; title('AP107 mPFC post-learning ROIs');

subplot(2,2,1)
plot(timevec, original_left_frontal_ROIs(1,:)', 'r'); 
xline(0.15)
ylabel('Left ROI fluorescence')
% xlabel('Time from stim onset')
ylim([-3*10^-3 7*10^-3])

subplot(2,2,3)
plot(timevec, original_right_frontal_ROIs(1,:)', 'blue'); 
xline(0.15)
ylabel('Right ROI fluorescence')
% xlabel('Time from stim onset')
ylim([-3*10^-3 7*10^-3])

subplot(2,2,2)
plot(timevec, reversal_left_frontal_ROIs(1,:)', 'r'); 
xline(0.15)
% ylabel('Left ROI fluorescence')
xlabel('Time from stim onset')
ylim([-3*10^-3 7*10^-3])

subplot(2,2,4)
plot(timevec, reversal_right_frontal_ROIs(1,:)', 'blue'); 
xline(0.15)
% ylabel('Right ROI fluorescence')
xlabel('Time from stim onset')
ylim([-3*10^-3 7*10^-3])

% all together
roi_colors = cat(1,brewermap(2,'Reds'),brewermap(2,'Blues'));

figure;
title('mPFC ROIs in original and reversal task learner')
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

%% Scroll movie for pre-learning (first:day) and post-learning (day:end)

% load activity from task
load('task_master_activity.mat');
animals = {task_master_activity.animal};

timestep = task_master_activity.timestep;
start_time = task_master_activity.start_time;
end_time = task_master_activity.end_time;
stim_frame = task_master_activity.stim_frame;
timevec = task_master_activity.timevec;
num_comp = task_master_activity.num_comp;

% choose averaging window
small_window = [0.05 0.3];
small_window_idx = timevec>=small_window(1)&timevec<=small_window(2);

% load master U
master_u_fn = 'D:\Andrada\Master project\widefield_alignment\U_master';
load(master_u_fn);

% possible stims
possible_stimuli = [-1 0 1];

%% - Original task

% define learned day
learned_day_idx = 5;

%% -- pre-learning
all_mice_avg_act = nan(num_comp,length(timevec),length(animals));

for animal_id = 1:length(animals)
    % get day indexes
    all_training_days = task_master_activity(animal_id).day;
    original_task_days_mask = task_master_activity(animal_id).original_task_days_mask;
    muscimol_days_mask = task_master_activity(animal_id).muscimol_days_mask;
    original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);
    
    % find which day the learned day is in whole dataset
    these_day_idx = original_task_day_idx(1:learned_day_idx-1);
    
    this_mouse_act = nan(num_comp, length(timevec), length(these_day_idx));
    for day_idx = these_day_idx
        % get stim act for this animal and day
        stim_act = task_master_activity(animal_id).stim_activity{day_idx};
        
        % load trial information and wheel move
        trialStimulusValue = task_master_activity(animal_id).trial_id{day_idx};
        % only completed trials
        trialStimulusValue = trialStimulusValue(1:size(stim_act,3));
        
        stim_wheel_move = task_master_activity(animal_id).stim_wheel_move{day_idx};
        
        % find average activity for each stimulus
        stim_avg_act = nanmean(stim_act,3);
        stim_avg_act = stim_avg_act - stim_avg_act(:,stim_frame,:);
        
        % deconvolve
        deconvolved_stim_avg_act = AP_deconv_wf(stim_avg_act, [], 1/timestep);
        
        this_mouse_act(:,:,day_idx-these_day_idx(1)+1) = deconvolved_stim_avg_act;
    end
    
    all_mice_avg_act(:,:,animal_id) = mean(this_mouse_act,3);
end

% get fluorescence
pre_learned_avg_act = mean(all_mice_avg_act,3);
pre_learned_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),pre_learned_avg_act);

% video
AP_image_scroll(pre_learned_avg_fluorescence,timevec); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)])
axis image;

%% -- post-learning
all_mice_avg_act = nan(num_comp,length(timevec),length(animals));

for animal_id = 1:length(animals)
    % get day indexes
    all_training_days = task_master_activity(animal_id).day;
    original_task_days_mask = task_master_activity(animal_id).original_task_days_mask;
    muscimol_days_mask = task_master_activity(animal_id).muscimol_days_mask;
    original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);
    
    % find which day the learned day is in whole dataset
    these_day_idx = original_task_day_idx(learned_day_idx:end);
    
    this_mouse_act = nan(num_comp, length(timevec), length(these_day_idx));
    for day_idx = these_day_idx
        % get stim act for this animal and day
        stim_act = task_master_activity(animal_id).stim_activity{day_idx};
        
        % load trial information and wheel move
        trialStimulusValue = task_master_activity(animal_id).trial_id{day_idx};
        % only completed trials
        trialStimulusValue = trialStimulusValue(1:size(stim_act,3));
        
        stim_wheel_move = task_master_activity(animal_id).stim_wheel_move{day_idx};
        
        % find average activity for each stimulus
        stim_avg_act = nanmean(stim_act,3);
        stim_avg_act = stim_avg_act - stim_avg_act(:,stim_frame,:);
        
        % deconvolve
        deconvolved_stim_avg_act = AP_deconv_wf(stim_avg_act, [], 1/timestep);
        
        this_mouse_act(:,:,day_idx-these_day_idx(1)+1) = deconvolved_stim_avg_act;
    end
    
    all_mice_avg_act(:,:,animal_id) = mean(this_mouse_act,3);
end

% get fluorescence
post_learned_avg_act = mean(all_mice_avg_act,3);
post_learned_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),post_learned_avg_act);

% video
AP_image_scroll(post_learned_avg_fluorescence,timevec); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)])
axis image;

%% - Reversal task

% define learned day
learned_day_idx = 7;

%% - non-learners
%% -- pre-learning
all_mice_avg_act = nan(num_comp,length(timevec),length(animals));

for animal_id = 2:length(animals)
    % get day indexes
    all_training_days = task_master_activity(animal_id).day;
    reversal_task_days_mask = task_master_activity(animal_id).reversal_task_days_mask;
    reversal_task_day_idx = find(reversal_task_days_mask);
    
    % find which day the learned day is in whole dataset
    these_day_idx = reversal_task_day_idx(1:learned_day_idx-1);
    
    this_mouse_act = nan(num_comp, length(timevec), length(these_day_idx));
    for day_idx = these_day_idx
        % get stim act for this animal and day
        stim_act = task_master_activity(animal_id).stim_activity{day_idx};
        
        % load trial information and wheel move
        trialStimulusValue = task_master_activity(animal_id).trial_id{day_idx};
        % only completed trials
        trialStimulusValue = trialStimulusValue(1:size(stim_act,3));
        
        stim_wheel_move = task_master_activity(animal_id).stim_wheel_move{day_idx};
        
        % find average activity for each stimulus
        stim_avg_act = nanmean(stim_act,3);
        stim_avg_act = stim_avg_act - stim_avg_act(:,stim_frame,:);
        
        % deconvolve
        deconvolved_stim_avg_act = AP_deconv_wf(stim_avg_act, [], 1/timestep);
        
        this_mouse_act(:,:,day_idx-these_day_idx(1)+1) = deconvolved_stim_avg_act;
    end
    
    all_mice_avg_act(:,:,animal_id) = mean(this_mouse_act,3);
end

% get fluorescence
pre_learned_avg_act = nanmean(all_mice_avg_act,3);
pre_learned_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),pre_learned_avg_act);

% video
AP_image_scroll(pre_learned_avg_fluorescence,timevec); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)])
axis image;

%% -- post-learning
all_mice_avg_act = nan(num_comp,length(timevec),length(animals));

for animal_id = 2:length(animals)
    % get day indexes
    all_training_days = task_master_activity(animal_id).day;
    reversal_task_days_mask = task_master_activity(animal_id).reversal_task_days_mask;
    reversal_task_day_idx = find(reversal_task_days_mask);
    
    % find which day the learned day is in whole dataset
    these_day_idx = reversal_task_day_idx(learned_day_idx:end);
    
    this_mouse_act = nan(num_comp, length(timevec), length(these_day_idx));
    for day_idx = these_day_idx
        % get stim act for this animal and day
        stim_act = task_master_activity(animal_id).stim_activity{day_idx};
        
        % load trial information and wheel move
        trialStimulusValue = task_master_activity(animal_id).trial_id{day_idx};
        % only completed trials
        trialStimulusValue = trialStimulusValue(1:size(stim_act,3));
        
        stim_wheel_move = task_master_activity(animal_id).stim_wheel_move{day_idx};
        
        % find average activity for each stimulus
        stim_avg_act = nanmean(stim_act,3);
        stim_avg_act = stim_avg_act - stim_avg_act(:,stim_frame,:);
        
        % deconvolve
        deconvolved_stim_avg_act = AP_deconv_wf(stim_avg_act, [], 1/timestep);
        
        this_mouse_act(:,:,day_idx-these_day_idx(1)+1) = deconvolved_stim_avg_act;
    end
    
    all_mice_avg_act(:,:,animal_id) = mean(this_mouse_act,3);
end

% get fluorescence
post_learned_avg_act = nanmean(all_mice_avg_act,3);
post_learned_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),post_learned_avg_act);

% video
AP_image_scroll(post_learned_avg_fluorescence,timevec); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)])
axis image;

%% learner
%% -- pre-learning

animal_id = 1;

% get day indexes
all_training_days = task_master_activity(animal_id).day;
reversal_task_days_mask = task_master_activity(animal_id).reversal_task_days_mask;
reversal_task_day_idx = find(reversal_task_days_mask);

% find which day the learned day is in whole dataset
these_day_idx = reversal_task_day_idx(1:learned_day_idx-1);

this_mouse_act = nan(num_comp, length(timevec), length(these_day_idx));
for day_idx = these_day_idx
    % get stim act for this animal and day
    stim_act = task_master_activity(animal_id).stim_activity{day_idx};
    
    % load trial information and wheel move
    trialStimulusValue = task_master_activity(animal_id).trial_id{day_idx};
    % only completed trials
    trialStimulusValue = trialStimulusValue(1:size(stim_act,3));
    
    stim_wheel_move = task_master_activity(animal_id).stim_wheel_move{day_idx};
    
    % find average activity for each stimulus
    stim_avg_act = nanmean(stim_act,3);
    stim_avg_act = stim_avg_act - stim_avg_act(:,stim_frame,:);
    
    % deconvolve
    deconvolved_stim_avg_act = AP_deconv_wf(stim_avg_act, [], 1/timestep);
    
    this_mouse_act(:,:,day_idx-these_day_idx(1)+1) = deconvolved_stim_avg_act;
end

all_mice_avg_act = mean(this_mouse_act,3);

% get fluorescence
pre_learned_avg_act = mean(all_mice_avg_act,3);
pre_learned_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),pre_learned_avg_act);

% video
AP_image_scroll(pre_learned_avg_fluorescence,timevec); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)])
axis image;

%% -- post-learning

% get day indexes
all_training_days = task_master_activity(animal_id).day;
reversal_task_days_mask = task_master_activity(animal_id).reversal_task_days_mask;
reversal_task_day_idx = find(reversal_task_days_mask);

% find which day the learned day is in whole dataset
these_day_idx = reversal_task_day_idx(learned_day_idx:end);

this_mouse_act = nan(num_comp, length(timevec), length(these_day_idx));
for day_idx = these_day_idx
    % get stim act for this animal and day
    stim_act = task_master_activity(animal_id).stim_activity{day_idx};
    
    % load trial information and wheel move
    trialStimulusValue = task_master_activity(animal_id).trial_id{day_idx};
    % only completed trials
    trialStimulusValue = trialStimulusValue(1:size(stim_act,3));
        
    % find average activity for each stimulus
    stim_avg_act = nanmean(stim_act,3);
    stim_avg_act = stim_avg_act - stim_avg_act(:,stim_frame,:);
    
    % deconvolve
    deconvolved_stim_avg_act = AP_deconv_wf(stim_avg_act, [], 1/timestep);
    
    this_mouse_act(:,:,day_idx-these_day_idx(1)+1) = deconvolved_stim_avg_act;
end

all_mice_avg_act = mean(this_mouse_act,3);

% get fluorescence
post_learned_avg_act = mean(all_mice_avg_act,3);
post_learned_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),post_learned_avg_act);

% video
AP_image_scroll(post_learned_avg_fluorescence,timevec); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)])
axis image;

%% - save movies for learner pre and post
% pre
AP_image_scroll(pre_learned_avg_fluorescence,timevec); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)])
axis image;
im =pre_learned_avg_fluorescence;
framerate = 50;
color_map = get(gcf,'colormap');
color_axis = caxis;
figure_position = get(gcf, 'position');
savefile = 'pre_learn_movie.avi';
t_annotation = timevec;
AP_movie2avi(im,framerate,color_map,color_axis,figure_position,savefile)

% post
AP_image_scroll(post_learned_avg_fluorescence,timevec); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)])
axis image;
im =post_learned_avg_fluorescence;
framerate = 50;
color_map = get(gcf,'colormap');
color_axis = caxis;
figure_position = get(gcf, 'position');
savefile = 'post_learn_movie.avi';
t_annotation = timevec;
AP_movie2avi(im,framerate,color_map,color_axis,figure_position,savefile)


%% Diff between original and reversal post-learning 

%% -- post-learning original

all_mice_avg_act = nan(num_comp,length(timevec),length(animals));

for animal_id = 1:length(animals)
    % get day indexes
    all_training_days = task_master_activity(animal_id).day;
    original_task_days_mask = task_master_activity(animal_id).original_task_days_mask;
    muscimol_days_mask = task_master_activity(animal_id).muscimol_days_mask;
    original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);
    
    % find which day the learned day is in whole dataset
    these_day_idx = original_task_day_idx(learned_day_idx:end);
    
    this_mouse_act = nan(num_comp, length(timevec), length(these_day_idx));
    for day_idx = these_day_idx
        % get stim act for this animal and day
        stim_act = task_master_activity(animal_id).stim_activity{day_idx};
        
        % load trial information and wheel move
        trialStimulusValue = task_master_activity(animal_id).trial_id{day_idx};
        % only completed trials
        trialStimulusValue = trialStimulusValue(1:size(stim_act,3));
        
        stim_wheel_move = task_master_activity(animal_id).stim_wheel_move{day_idx};
        
        % find average activity for each stimulus
        stim_avg_act = nanmean(stim_act,3);
        stim_avg_act = stim_avg_act - stim_avg_act(:,stim_frame,:);
        
        % deconvolve
        deconvolved_stim_avg_act = AP_deconv_wf(stim_avg_act, [], 1/timestep);
        
        this_mouse_act(:,:,day_idx-these_day_idx(1)+1) = deconvolved_stim_avg_act;
    end
    
    all_mice_avg_act(:,:,animal_id) = mean(this_mouse_act,3);
end

% get 107 data
original_avg_act = all_mice_avg_act(:,:,1);
original_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),original_avg_act);

% get for others
other_original_avg_act = mean(all_mice_avg_act(:,:,2:end),3);
other_original_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),other_original_avg_act);

% video
AP_image_scroll(original_avg_fluorescence,timevec); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)])
axis image;

AP_image_scroll(other_original_avg_fluorescence,timevec); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)])
axis image;


%% - post-learning reversal

all_mice_avg_act = nan(num_comp,length(timevec),length(animals));

for animal_id = 1:length(animals)
    % get day indexes
    all_training_days = task_master_activity(animal_id).day;
    reversal_task_days_mask = task_master_activity(animal_id).reversal_task_days_mask;
    reversal_task_day_idx = find(reversal_task_days_mask);
    
    % find which day the learned day is in whole dataset
    these_day_idx = reversal_task_day_idx(learned_day_idx:end);
    
    this_mouse_act = nan(num_comp, length(timevec), length(these_day_idx));
    for day_idx = these_day_idx
        % get stim act for this animal and day
        stim_act = task_master_activity(animal_id).stim_activity{day_idx};
        
        % load trial information and wheel move
        trialStimulusValue = task_master_activity(animal_id).trial_id{day_idx};
        % only completed trials
        trialStimulusValue = trialStimulusValue(1:size(stim_act,3));
        
        stim_wheel_move = task_master_activity(animal_id).stim_wheel_move{day_idx};
        
        % find average activity for each stimulus
        stim_avg_act = nanmean(stim_act,3);
        stim_avg_act = stim_avg_act - stim_avg_act(:,stim_frame,:);
        
        % deconvolve
        deconvolved_stim_avg_act = AP_deconv_wf(stim_avg_act, [], 1/timestep);
        
        this_mouse_act(:,:,day_idx-these_day_idx(1)+1) = deconvolved_stim_avg_act;
    end
    
    all_mice_avg_act(:,:,animal_id) = mean(this_mouse_act,3);
end

% get fluorescence for 107
reversal_avg_act = all_mice_avg_act(:,:,1);
reversal_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),reversal_avg_act);

% get for others
other_reversal_avg_act = mean(all_mice_avg_act(:,:,2:end),3);
other_reversal_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),other_reversal_avg_act);

% videos
AP_image_scroll(reversal_avg_fluorescence,timevec); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)])
axis image;

AP_image_scroll(other_reversal_avg_fluorescence,timevec); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)])
axis image;

%% -- difference 
% -- difference for learner
diff_avg_act = reversal_avg_act - original_avg_act; 
diff_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),diff_avg_act);

% -- difference for non-learners
diff_other_avg_act = other_reversal_avg_act - other_original_avg_act; 
diff_other_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),diff_other_avg_act);

% -- difference between 
diff_between_avg_act = diff_avg_act - diff_other_avg_act;
diff_between_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),diff_between_avg_act);

% videos
AP_image_scroll(diff_avg_fluorescence,timevec); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)])
axis image;

AP_image_scroll(diff_other_avg_fluorescence,timevec); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)])
axis image;

AP_image_scroll(diff_between_avg_fluorescence,timevec); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)])
axis image;


% average brain in time window

% choose averaging window
small_window = [0.05 0.2];
small_window_idx = timevec>=small_window(1)&timevec<=small_window(2);

% average and get fluorescence
avg_diff_avg_act = mean(diff_avg_act(:,small_window_idx),2);
avg_diff_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),avg_diff_avg_act);

% -- difference for non-learners
avg_diff_other_avg_act = mean(diff_other_avg_act(:,small_window_idx),2);
avg_diff_other_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),avg_diff_other_avg_act);

% -- difference between 
avg_diff_between_avg_act = mean(diff_between_avg_act(:,small_window_idx),2);
avg_diff_between_avg_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),avg_diff_between_avg_act);

% plot of average brain for differences
figure; 
imagesc(avg_diff_avg_fluorescence); colormap(brewermap([], 'PRGn')); caxis([-5*10^(-3) 5*10^(-3)]); axis image;
title('Fluorescence difference reversal and original for learner (5-200 ms)');
set(gca,'Xtick', []);
set(gca,'Ytick', []);

figure; 
imagesc(avg_diff_other_avg_fluorescence); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)]); axis image;
title('Fluorescence difference reversal and original for non-learners (5-200 ms)');
set(gca,'Xtick', []);
set(gca,'Ytick', []);

figure; 
imagesc(avg_diff_between_avg_fluorescence); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)]); axis image;
title('Fluorescence difference reversal and original for learner minus difference for non-learners (5-200 ms)');
set(gca,'Xtick', []);
set(gca,'Ytick', []);

%% Brain plots avg 100ms after stim/move

small_window = [0 0.075];
small_window_idx = timevec>=small_window(1)&timevec<=small_window(2);

%% - Stim Original task

% define learned day
learned_day_idx = 5;

%% -- pre-learning
all_mice_avg_act = nan(num_comp,length(timevec),length(animals));

for animal_id = 1:length(animals)
    % get day indexes
    all_training_days = task_master_activity(animal_id).day;
    original_task_days_mask = task_master_activity(animal_id).original_task_days_mask;
    muscimol_days_mask = task_master_activity(animal_id).muscimol_days_mask;
    original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);
    
    % find which day the learned day is in whole dataset
    these_day_idx = original_task_day_idx(1:learned_day_idx-1);
    
    this_mouse_act = nan(num_comp, length(timevec), length(these_day_idx));
    for day_idx = these_day_idx
        % get stim act for this animal and day
        stim_act = task_master_activity(animal_id).stim_activity{day_idx};
        
        % load trial information and wheel move
        trialStimulusValue = task_master_activity(animal_id).trial_id{day_idx};
        % only completed trials
        trialStimulusValue = trialStimulusValue(1:size(stim_act,3));
                
        % find average activity for each stimulus
        stim_avg_act = nanmean(stim_act,3);
        stim_avg_act = stim_avg_act - stim_avg_act(:,stim_frame,:);
        
        % deconvolve
        deconvolved_stim_avg_act = AP_deconv_wf(stim_avg_act, [], 1/timestep);
        
        this_mouse_act(:,:,day_idx-these_day_idx(1)+1) = deconvolved_stim_avg_act;
    end
    
    all_mice_avg_act(:,:,animal_id) = mean(this_mouse_act,3);
end

pre_learning_original_act = mean(all_mice_avg_act(:,small_window_idx,:),2);
pre_learning_original_act = permute(pre_learning_original_act, [1,3,2]);
pre_learning_original_act = mean(pre_learning_original_act, 2);
pre_learned_original_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),pre_learning_original_act);

%% -- post-learning
all_mice_avg_act = nan(num_comp,length(timevec),length(animals));

for animal_id = 1:length(animals)
    % get day indexes
    all_training_days = task_master_activity(animal_id).day;
    original_task_days_mask = task_master_activity(animal_id).original_task_days_mask;
    muscimol_days_mask = task_master_activity(animal_id).muscimol_days_mask;
    original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);
    
    % find which day the learned day is in whole dataset
    these_day_idx = original_task_day_idx(learned_day_idx:end);
    
    this_mouse_act = nan(num_comp, length(timevec), length(these_day_idx));
    for day_idx = these_day_idx
        % get stim act for this animal and day
        stim_act = task_master_activity(animal_id).stim_activity{day_idx};
        
        % load trial information and wheel move
        trialStimulusValue = task_master_activity(animal_id).trial_id{day_idx};
        % only completed trials
        trialStimulusValue = trialStimulusValue(1:size(stim_act,3));
        
        % find average activity for each stimulus
        stim_avg_act = nanmean(stim_act,3);
        stim_avg_act = stim_avg_act - stim_avg_act(:,stim_frame,:);
        
        % deconvolve
        deconvolved_stim_avg_act = AP_deconv_wf(stim_avg_act, [], 1/timestep);
        
        this_mouse_act(:,:,day_idx-these_day_idx(1)+1) = deconvolved_stim_avg_act;
    end
    
    all_mice_avg_act(:,:,animal_id) = mean(this_mouse_act,3);
end

post_learning_original_act = mean(all_mice_avg_act(:,small_window_idx,:),2);
post_learning_original_act = permute(post_learning_original_act, [1,3,2]);
post_learning_original_act = mean(post_learning_original_act, 2);
post_learned_original_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),post_learning_original_act);

%% - Stim Reversal task

% define learned day
learned_day_idx = 7;

%% -- pre-learning
all_mice_avg_act = nan(num_comp,length(timevec),length(animals));

for animal_id = 1:length(animals)
    % get day indexes
    all_training_days = task_master_activity(animal_id).day;
    reversal_task_days_mask = task_master_activity(animal_id).reversal_task_days_mask;
    reversal_task_day_idx = find(reversal_task_days_mask);
    
    % find which day the learned day is in whole dataset
    these_day_idx = reversal_task_day_idx(1:learned_day_idx-1);
    
    this_mouse_act = nan(num_comp, length(timevec), length(these_day_idx));
    for day_idx = these_day_idx
        % get stim act for this animal and day
        stim_act = task_master_activity(animal_id).stim_activity{day_idx};
        
        % load trial information and wheel move
        trialStimulusValue = task_master_activity(animal_id).trial_id{day_idx};
        % only completed trials
        trialStimulusValue = trialStimulusValue(1:size(stim_act,3));
        
        stim_wheel_move = task_master_activity(animal_id).stim_wheel_move{day_idx};
        
        % find average activity for each stimulus
        stim_avg_act = nanmean(stim_act,3);
        stim_avg_act = stim_avg_act - stim_avg_act(:,stim_frame,:);
        
        % deconvolve
        deconvolved_stim_avg_act = AP_deconv_wf(stim_avg_act, [], 1/timestep);
        
        this_mouse_act(:,:,day_idx-these_day_idx(1)+1) = deconvolved_stim_avg_act;
    end
    
    all_mice_avg_act(:,:,animal_id) = mean(this_mouse_act,3);
end

pre_learning_reversal_act = mean(all_mice_avg_act(:,small_window_idx,1),2);
pre_learning_reversal_act = permute(pre_learning_reversal_act, [1,3,2]);
pre_learning_reversal_act = mean(pre_learning_reversal_act, 2);
pre_learned_reversal_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),pre_learning_reversal_act);

pre_learning_others_reversal_act = mean(all_mice_avg_act(:,small_window_idx,2:end),2);
pre_learning_others_reversal_act = permute(pre_learning_others_reversal_act, [1,3,2]);
pre_learning_others_reversal_act = mean(pre_learning_others_reversal_act, 2);
pre_learned_others_reversal_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),pre_learning_others_reversal_act);

%% -- post-learning
all_mice_avg_act = nan(num_comp,length(timevec),length(animals));

for animal_id = 1:length(animals)
    % get day indexes
    all_training_days = task_master_activity(animal_id).day;
    reversal_task_days_mask = task_master_activity(animal_id).reversal_task_days_mask;
    reversal_task_day_idx = find(reversal_task_days_mask);
    
    % find which day the learned day is in whole dataset
    these_day_idx = reversal_task_day_idx(learned_day_idx:end);
    
    this_mouse_act = nan(num_comp, length(timevec), length(these_day_idx));
    for day_idx = these_day_idx
        % get stim act for this animal and day
        stim_act = task_master_activity(animal_id).stim_activity{day_idx};
        
        % load trial information and wheel move
        trialStimulusValue = task_master_activity(animal_id).trial_id{day_idx};
        % only completed trials
        trialStimulusValue = trialStimulusValue(1:size(stim_act,3));
        
        stim_wheel_move = task_master_activity(animal_id).stim_wheel_move{day_idx};
        
        % find average activity for each stimulus
        stim_avg_act = nanmean(stim_act,3);
        stim_avg_act = stim_avg_act - stim_avg_act(:,stim_frame,:);
        
        % deconvolve
        deconvolved_stim_avg_act = AP_deconv_wf(stim_avg_act, [], 1/timestep);
        
        this_mouse_act(:,:,day_idx-these_day_idx(1)+1) = deconvolved_stim_avg_act;
    end
    
    all_mice_avg_act(:,:,animal_id) = mean(this_mouse_act,3);
end

post_learning_reversal_act = mean(all_mice_avg_act(:,small_window_idx,1),2);
post_learning_reversal_act = permute(post_learning_reversal_act, [1,3,2]);
post_learning_reversal_act = mean(post_learning_reversal_act, 2);
post_learned_reversal_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),post_learning_reversal_act);

post_learning_others_reversal_act = mean(all_mice_avg_act(:,small_window_idx,2:end),2);
post_learning_others_reversal_act = permute(post_learning_others_reversal_act, [1,3,2]);
post_learning_others_reversal_act = mean(post_learning_others_reversal_act, 2);
post_learned_others_reversal_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),post_learning_others_reversal_act);

%% - Stim plot

figure;
subplot(2,3,1)
imagesc(pre_learned_original_fluorescence); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)]); axis image;
title('Original task');
set(gca,'Xtick', []);
set(gca,'Ytick', []);

subplot(2,3,4)
imagesc(post_learned_original_fluorescence); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)]); axis image;
set(gca,'Xtick', []);
set(gca,'Ytick', []);

subplot(2,3,2)
imagesc(pre_learned_reversal_fluorescence); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)]); axis image;
title('Reversal task learner');
set(gca,'Xtick', []);
set(gca,'Ytick', []);

subplot(2,3,5)
imagesc(post_learned_reversal_fluorescence); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)]); axis image;
set(gca,'Xtick', []);
set(gca,'Ytick', []);

subplot(2,3,3)
imagesc(pre_learned_others_reversal_fluorescence); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)]); axis image;
title('Reversal task non-learners');
set(gca,'Xtick', []);
set(gca,'Ytick', []);

subplot(2,3,6)
imagesc(post_learned_others_reversal_fluorescence); colormap(brewermap([], 'PRGn')); caxis([-7*10^(-3) 7*10^(-3)]); axis image;
set(gca,'Xtick', []);
set(gca,'Ytick', []);
%% - Move Original task

small_window = [0 0.075];
small_window_idx = timevec>=small_window(1)&timevec<=small_window(2);

move_frame = task_master_activity.stim_frame;

% define learned day
learned_day_idx = 5;

%% -- pre-learning
all_mice_avg_act = nan(num_comp,length(timevec),length(animals));

for animal_id = 1:length(animals)
    % get day indexes
    all_training_days = task_master_activity(animal_id).day;
    original_task_days_mask = task_master_activity(animal_id).original_task_days_mask;
    muscimol_days_mask = task_master_activity(animal_id).muscimol_days_mask;
    original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);
    
    % find which day the learned day is in whole dataset
    these_day_idx = original_task_day_idx(1:learned_day_idx-1);
    
    this_mouse_act = nan(num_comp, length(timevec), length(these_day_idx));
    for day_idx = these_day_idx
        % get stim act for this animal and day
        move_act = task_master_activity(animal_id).move_activity{day_idx};
        
        % load trial information and wheel move
        trialStimulusValue = task_master_activity(animal_id).trial_id{day_idx};
        % only completed trials
        trialStimulusValue = trialStimulusValue(1:size(move_act,3));
                
        % find average activity for each stimulus
        move_avg_act = nanmean(move_act,3);
        move_avg_act = move_avg_act - move_avg_act(:,move_frame,:);
        
        % deconvolve
        deconvolved_move_avg_act = AP_deconv_wf(move_avg_act, [], 1/timestep);
        
        this_mouse_act(:,:,day_idx-these_day_idx(1)+1) = deconvolved_move_avg_act;
    end
    
    all_mice_avg_act(:,:,animal_id) = mean(this_mouse_act,3);
end

pre_learning_original_act = mean(all_mice_avg_act(:,small_window_idx,:),2);
pre_learning_original_act = permute(pre_learning_original_act, [1,3,2]);
pre_learning_original_act = mean(pre_learning_original_act, 2);
pre_learned_original_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),pre_learning_original_act);

%% -- post-learning
all_mice_avg_act = nan(num_comp,length(timevec),length(animals));

for animal_id = 1:length(animals)
    % get day indexes
    all_training_days = task_master_activity(animal_id).day;
    original_task_days_mask = task_master_activity(animal_id).original_task_days_mask;
    muscimol_days_mask = task_master_activity(animal_id).muscimol_days_mask;
    original_task_day_idx = find(original_task_days_mask-muscimol_days_mask==1);
    
    % find which day the learned day is in whole dataset
    these_day_idx = original_task_day_idx(learned_day_idx:end);
    
    this_mouse_act = nan(num_comp, length(timevec), length(these_day_idx));
    for day_idx = these_day_idx
        % get stim act for this animal and day
        move_act = task_master_activity(animal_id).move_activity{day_idx};
        
        % load trial information and wheel move
        trialStimulusValue = task_master_activity(animal_id).trial_id{day_idx};
        % only completed trials
        trialStimulusValue = trialStimulusValue(1:size(move_act,3));
        
        % find average activity for each stimulus
        move_avg_act = nanmean(move_act,3);
        move_avg_act = move_avg_act - move_avg_act(:,move_frame,:);
        
        % deconvolve
        deconvolved_move_avg_act = AP_deconv_wf(move_avg_act, [], 1/timestep);
        
        this_mouse_act(:,:,day_idx-these_day_idx(1)+1) = deconvolved_move_avg_act;
    end
    
    all_mice_avg_act(:,:,animal_id) = mean(this_mouse_act,3);
end

post_learning_original_act = mean(all_mice_avg_act(:,small_window_idx,:),2);
post_learning_original_act = permute(post_learning_original_act, [1,3,2]);
post_learning_original_act = mean(post_learning_original_act, 2);
post_learned_original_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),post_learning_original_act);

%% - Move Reversal task

% define learned day
learned_day_idx = 7;

%% -- pre-learning
all_mice_avg_act = nan(num_comp,length(timevec),length(animals));

for animal_id = 1:length(animals)
    % get day indexes
    all_training_days = task_master_activity(animal_id).day;
    reversal_task_days_mask = task_master_activity(animal_id).reversal_task_days_mask;
    reversal_task_day_idx = find(reversal_task_days_mask);
    
    % find which day the learned day is in whole dataset
    these_day_idx = reversal_task_day_idx(1:learned_day_idx-1);
    
    this_mouse_act = nan(num_comp, length(timevec), length(these_day_idx));
    for day_idx = these_day_idx
        % get stim act for this animal and day
        move_act = task_master_activity(animal_id).move_activity{day_idx};
        
        % load trial information and wheel move
        trialStimulusValue = task_master_activity(animal_id).trial_id{day_idx};
        % only completed trials
        trialStimulusValue = trialStimulusValue(1:size(move_act,3));
       
        % find average activity for each stimulus
        move_avg_act = nanmean(move_act,3);
        move_avg_act = move_avg_act - move_avg_act(:,move_frame,:);
        
        % deconvolve
        deconvolved_move_avg_act = AP_deconv_wf(move_avg_act, [], 1/timestep);
        
        this_mouse_act(:,:,day_idx-these_day_idx(1)+1) = deconvolved_move_avg_act;
    end
    
    all_mice_avg_act(:,:,animal_id) = mean(this_mouse_act,3);
end

pre_learning_reversal_act = mean(all_mice_avg_act(:,small_window_idx,1),2);
pre_learning_reversal_act = permute(pre_learning_reversal_act, [1,3,2]);
pre_learning_reversal_act = mean(pre_learning_reversal_act, 2);
pre_learned_reversal_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),pre_learning_reversal_act);

pre_learning_others_reversal_act = mean(all_mice_avg_act(:,small_window_idx,2:end),2);
pre_learning_others_reversal_act = permute(pre_learning_others_reversal_act, [1,3,2]);
pre_learning_others_reversal_act = mean(pre_learning_others_reversal_act, 2);
pre_learned_others_reversal_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),pre_learning_others_reversal_act);

%% -- post-learning
all_mice_avg_act = nan(num_comp,length(timevec),length(animals));

for animal_id = 1:length(animals)
    % get day indexes
    all_training_days = task_master_activity(animal_id).day;
    reversal_task_days_mask = task_master_activity(animal_id).reversal_task_days_mask;
    reversal_task_day_idx = find(reversal_task_days_mask);
    
    % find which day the learned day is in whole dataset
    these_day_idx = reversal_task_day_idx(learned_day_idx:end);
    
    this_mouse_act = nan(num_comp, length(timevec), length(these_day_idx));
    for day_idx = these_day_idx
        % get stim act for this animal and day
        move_act = task_master_activity(animal_id).move_activity{day_idx};
        
        % load trial information and wheel move
        trialStimulusValue = task_master_activity(animal_id).trial_id{day_idx};
        % only completed trials
        trialStimulusValue = trialStimulusValue(1:size(move_act,3));
        
        % find average activity for each stimulus
        move_avg_act = nanmean(move_act,3);
        move_avg_act = move_avg_act - move_avg_act(:,move_frame,:);
        
        % deconvolve
        deconvolved_move_avg_act = AP_deconv_wf(move_avg_act, [], 1/timestep);
        
        this_mouse_act(:,:,day_idx-these_day_idx(1)+1) = deconvolved_move_avg_act;
    end
    
    all_mice_avg_act(:,:,animal_id) = mean(this_mouse_act,3);
end

post_learning_reversal_act = mean(all_mice_avg_act(:,small_window_idx,1),2);
post_learning_reversal_act = permute(post_learning_reversal_act, [1,3,2]);
post_learning_reversal_act = mean(post_learning_reversal_act, 2);
post_learned_reversal_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),post_learning_reversal_act);

post_learning_others_reversal_act = mean(all_mice_avg_act(:,small_window_idx,2:end),2);
post_learning_others_reversal_act = permute(post_learning_others_reversal_act, [1,3,2]);
post_learning_others_reversal_act = mean(post_learning_others_reversal_act, 2);
post_learned_others_reversal_fluorescence = AP_svdFrameReconstruct(U_master(:,:,1:num_comp),post_learning_others_reversal_act);

%% - Move plot

figure('Name', 'Move plot');
subplot(2,3,1)
imagesc(pre_learned_original_fluorescence); colormap(brewermap([], 'PRGn')); caxis([-14*10^(-3) 14*10^(-3)]); axis image;
title('Original task');
set(gca,'Xtick', []);
set(gca,'Ytick', []);

subplot(2,3,4)
imagesc(post_learned_original_fluorescence); colormap(brewermap([], 'PRGn')); caxis([-14*10^(-3) 14*10^(-3)]); axis image;
set(gca,'Xtick', []);
set(gca,'Ytick', []);

subplot(2,3,2)
imagesc(pre_learned_reversal_fluorescence); colormap(brewermap([], 'PRGn')); caxis([-14*10^(-3) 14*10^(-3)]); axis image;
title('Reversal task learner');
set(gca,'Xtick', []);
set(gca,'Ytick', []);

subplot(2,3,5)
imagesc(post_learned_reversal_fluorescence); colormap(brewermap([], 'PRGn')); caxis([-14*10^(-3) 14*10^(-3)]); axis image;
set(gca,'Xtick', []);
set(gca,'Ytick', []);

subplot(2,3,3)
imagesc(pre_learned_others_reversal_fluorescence); colormap(brewermap([], 'PRGn')); caxis([-14*10^(-3) 14*10^(-3)]); axis image;
title('Reversal task non-learners');
set(gca,'Xtick', []);
set(gca,'Ytick', []);

subplot(2,3,6)
imagesc(post_learned_others_reversal_fluorescence); colormap(brewermap([], 'PRGn')); caxis([-14*10^(-3) 14*10^(-3)]); axis image;
set(gca,'Xtick', []);
set(gca,'Ytick', []);