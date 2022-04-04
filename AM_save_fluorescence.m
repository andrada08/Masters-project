%% Get and save passive fluorescence

animals = {'AP107','AP108', 'AP113', 'AP114', 'AP115'};

master_u_fn = 'D:\Andrada\Master project\widefield_alignment\U_master';
load(master_u_fn);

passive_master_activity = struct;

% create matrix of times for movie
timestep = 0.01;
start_time = -0.5;
end_time = 1;
timevec = [start_time:timestep:end_time];
stim_frame = (-start_time)*(1/timestep)+1;


passive_master_activity.timestep = timestep;
passive_master_activity.start_time = start_time;
passive_master_activity.end_time = end_time;
passive_master_activity.stim_frame = stim_frame;
passive_master_activity.timevec = timevec;

% define number of components
num_comp = 200;
passive_master_activity.num_comp = 200;

for animal_id=1:length(animals)
    animal = animals{animal_id};
    passive_master_activity(animal_id).animal = animal;
    
    % find passive days
    protocol = 'AP_lcrGratingPassive';
    experiments = AP_find_experiments(animal,protocol);
    for day_index=1:length({experiments.day})
        day = experiments(day_index).day;
        passive_master_activity(animal_id).day{day_index} = day;
        
        % find passive experiment
        experiment = experiments(day_index).experiment(end);
        
        % load experiment
        load_parts.imaging = true;
        load_parts.cam = true;
        verbose = true;
        AP_load_experiment;
        
        % find value of stimulus per trial
        trialStimulusValue = signals_events.stimAzimuthValues/90 .* signals_events.stimContrastValues;
        
        % align U's
        Udf_aligned = AP_align_widefield(Udf,animal,day);
        % get corresponding V's
        fVdf_Umaster = ChangeU(Udf_aligned,fVdf,U_master);
                       
        % find activity for the time defined
        time_stimulus = stimOn_times+timevec;
        all_stim_act = interp1(frame_t,fVdf_Umaster',time_stimulus);
        all_stim_act = permute(all_stim_act, [3,2,1]);
        all_stim_act = all_stim_act - all_stim_act(:,stim_frame,:);
        
               
        % save in struct
        passive_master_activity(animal_id).stim_activity{day_index} = all_stim_act(1:num_comp,:,:);
        passive_master_activity(animal_id).trial_id{day_index} = trialStimulusValue;
        
        % stim aligned wheel move
        stim_wheel_move = interp1(Timeline.rawDAQTimestamps,+wheel_move,time_stimulus');
        passive_master_activity(animal_id).stim_wheel_move{day_index} = stim_wheel_move;

    end
    disp(['Done with ' animal])
end

disp('Done all')

% Save 
save('passive_master_activity.mat', 'passive_master_activity', '-v7.3')
disp('Saved')

%% Get and save task fluorescence

animals = {'AP107','AP108', 'AP113', 'AP114', 'AP115'};

master_u_fn = 'D:\Andrada\Master project\widefield_alignment\U_master';
load(master_u_fn);

task_master_activity = struct;

% create matrix of times for movie
timestep = 0.01;
start_time = -0.5;
end_time = 2;
timevec = [start_time:timestep:end_time];
stim_frame = (-start_time)*(1/timestep)+1;


task_master_activity.timestep = timestep;
task_master_activity.start_time = start_time;
task_master_activity.end_time = end_time;
task_master_activity.stim_frame = stim_frame;
task_master_activity.timevec = timevec;

% define number of components
num_comp = 200;
task_master_activity.num_comp = 200;

for animal_id=1:length(animals)
    animal = animals{animal_id};
    task_master_activity(animal_id).animal = animal;
    
    % find task days
    protocol = 'AP_stimWheel';
    flexible_name = true;
    experiments = AP_find_experiments(animal,protocol,flexible_name);
    for day_index=1:length({experiments.day})
        day = experiments(day_index).day;
        task_master_activity(animal_id).day{day_index} = day;
        
        % find task experiment
        experiment = experiments(day_index).experiment(end);
        
        % load experiment
        load_parts.imaging = true;
        load_parts.cam = true;
        verbose = true;
        AP_load_experiment;
        
        % find value of stimulus per trial
        trialStimulusValue = signals_events.trialContrastValues .* signals_events.trialSideValues;
        
        % align U's
        Udf_aligned = AP_align_widefield(Udf,animal,day);
        % get corresponding V's
        fVdf_Umaster = ChangeU(Udf_aligned,fVdf,U_master);
                       
        % find activity for the time defined
        time_stimulus = stimOn_times+timevec;
        all_stim_act = interp1(frame_t,fVdf_Umaster',time_stimulus);
        all_stim_act = permute(all_stim_act, [3,2,1]);
        all_stim_act = all_stim_act - all_stim_act(:,stim_frame,:);
        
               
        % save in struct
        task_master_activity(animal_id).stim_activity{day_index} = all_stim_act(1:num_comp,:,:);
        task_master_activity(animal_id).trial_id{day_index} = trialStimulusValue;
        
        % stim aligned wheel move
        stim_wheel_move = interp1(Timeline.rawDAQTimestamps,+wheel_move,time_stimulus');
        task_master_activity(animal_id).stim_wheel_move{day_index} = stim_wheel_move;

    end
    disp(['Done with ' animal])
end

disp('Done all')

% Save 
save('task_master_activity.mat', 'task_master_activity', '-v7.3')
disp('Saved')

%% Get all avg fluorescence and save in struct

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
end

save('passive_master_activity.mat', 'passive_master_activity', '-v7.3')
disp('Saved')

%% Get avg for each stim
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
    
    avg_frontal_left = nan(3,length(passive_master_activity(animal_id).day));
    avg_frontal_right = nan(3,length(passive_master_activity(animal_id).day));
    
    for stim_idx =1:length(possible_stimuli)
        for day_idx=1:length(passive_master_activity(animal_id).day)
            
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
            
            % activity for current stim
            this_stim_act = deconvolved_stim_avg_act(:,:,trialStimulusValue==possible_stimuli(stim_idx));
            
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
    end
end

save('passive_master_activity.mat', 'passive_master_activity', '-v7.3')
disp('Saved')

%% Get masks for passive only days, original task days and reversal task days

animals = {'AP107','AP108', 'AP113', 'AP114', 'AP115'};

load('passive_master_activity.mat');

for animal_id=1:length(animals)
    
    animal = animals{animal_id};
    
    % find muscimol days and save
    protocol = 'AP_sparseNoise';
    experiments = AP_find_experiments(animal,protocol);
    muscimol_days = {experiments(2:end).day};
    muscimol_days_mask = ismember(passive_master_activity(animal_id).day, muscimol_days);
    passive_master_activity(animal_id).muscimol_days_mask = muscimol_days_mask;
    
    % find indices for days of reversal task
    protocol = 'AP_stimWheelLeftReverse';
    experiments = AP_find_experiments(animal,protocol);
    reversal_task_days = {experiments.day};
    reversal_task_days_mask = ismember(passive_master_activity(animal_id).day, reversal_task_days);
    passive_master_activity(animal_id).reversal_task_days_mask = reversal_task_days_mask;
    
    % find indices for days for original task
    protocol = 'AP_stimWheelRight';
    experiments = AP_find_experiments(animal,protocol);
    original_task_days = {experiments.day};
    original_task_days_mask = ismember(passive_master_activity(animal_id).day, original_task_days);
    passive_master_activity(animal_id).original_task_days_mask = original_task_days_mask;
    
    % find indices for passive only days - no task and no muscimol
    passive_only_days_mask = (~ismember(passive_master_activity(animal_id).day, original_task_days))...
        &(~ismember(passive_master_activity(animal_id).day, reversal_task_days))...
        &(~ismember(passive_master_activity(animal_id).day, muscimol_days));
    passive_master_activity(animal_id).passive_only_days_mask = passive_only_days_mask;
end

save('passive_master_activity.mat', 'passive_master_activity', '-v7.3')
disp('Saved')

%% Get masks for original task days and reversal task days

animals = {'AP107','AP108', 'AP113', 'AP114', 'AP115'};

load('task_master_activity.mat');

for animal_id=1:length(animals)
    
    animal = animals{animal_id};
    
    % find muscimol days and save
    protocol = 'AP_sparseNoise';
    experiments = AP_find_experiments(animal,protocol);
    muscimol_days = {experiments(2:end).day};
    muscimol_days_mask = ismember(task_master_activity(animal_id).day, muscimol_days);
    task_master_activity(animal_id).muscimol_days_mask = muscimol_days_mask;
    
    % find indices for days of reversal task
    protocol = 'AP_stimWheelLeftReverse';
    experiments = AP_find_experiments(animal,protocol);
    reversal_task_days = {experiments.day};
    reversal_task_days_mask = ismember(task_master_activity(animal_id).day, reversal_task_days);
    task_master_activity(animal_id).reversal_task_days_mask = reversal_task_days_mask;
    
    % find indices for days for original task
    protocol = 'AP_stimWheelRight';
    experiments = AP_find_experiments(animal,protocol);
    original_task_days = {experiments.day};
    original_task_days_mask = ismember(task_master_activity(animal_id).day, original_task_days);
    task_master_activity(animal_id).original_task_days_mask = original_task_days_mask;
    
end

save('task_master_activity.mat', 'task_master_activity', '-v7.3')
disp('Saved')
