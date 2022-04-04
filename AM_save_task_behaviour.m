%% Save task behaviour for AP107

animals = {'AP107','AP108', 'AP113', 'AP114', 'AP115'};

behaviour = struct;

% create matrix of times 
timestep = 0.01;
start_time = -2;
end_time = 2;
timevec = [start_time:timestep:end_time];


behaviour.timestep = timestep;
behaviour.start_time = start_time;
behaviour.end_time = end_time;
behaviour.timevec = timevec;

for animal_id=1:length(animals)
    animal = animals{animal_id};
    behaviour(animal_id).animal = animal;
    
    % find muscimol days and save
    protocol = 'AP_sparseNoise';
    experiments = AP_find_experiments(animal,protocol);
    muscimol_days = {experiments(2:end).day};
    behaviour(animal_id).muscimol_days = muscimol_days;
    
    % find task days    
    protocol = 'AP_stimWheel';
    flexible_name = true;
    experiments = AP_find_experiments(animal,protocol,flexible_name);
    for day_index=1:length(experiments)
        day = experiments(day_index).day;
                
%         % skip if it's a muscimol day
%         is_muscimol = find(contains(muscimol_days,day));
%         if is_muscimol
%             continue
%         end

        % save if it's not
        behaviour(animal_id).day{day_index} = day;
        experiment = experiments(day_index).experiment(end);
        
        % load experiment
        load_parts.imaging = false;
        load_parts.cam = true;
        verbose = true;
        AP_load_experiment;
        
        % find value of stimulus per trial
        trialStimulusValue = signals_events.trialContrastValues .* signals_events.trialSideValues;
        behaviour(animal_id).trial_id{day_index} = trialStimulusValue;
        
        % time vec around stim onsets
        time_stimulus = stimOn_times+timevec;
        behaviour(animal_id).time_stimulus{day_index} = time_stimulus;
   
        % define t
        t = Timeline.rawDAQTimestamps;
        behaviour(animal_id).t{day_index} = t;
        
        % wheel position
        stim_wheel_position = interp1(t,wheel_position,time_stimulus');
        behaviour(animal_id).stim_wheel_position{day_index} = stim_wheel_position;
        
        % stim aligned wheel move
        stim_wheel_move = interp1(t,+wheel_move,time_stimulus');
        behaviour(animal_id).stim_wheel_move{day_index} = stim_wheel_move;
        
        % all moves
        tmp_move = [0; diff(wheel_move)];
        all_move_on_frames = find(tmp_move==1);
        all_move_on_times = t(all_move_on_frames);
        behaviour(animal_id).all_move_on_times{day_index} = all_move_on_times;
                
        % all move offsets 
        all_move_off_frames = find(tmp_move==-1);
        all_move_off_times = t(all_move_off_frames);
        behaviour(animal_id).all_move_off_times{day_index} = all_move_off_times;
        
        % move after stim times
        stim_move_on_times = all_move_on_times(cell2mat(arrayfun(@(X) find(all_move_on_times>X,1,'first'),stimOn_times', 'UniformOutput', 0)));
        behaviour(animal_id).stim_move_on_times{day_index} = stim_move_on_times;
        
        % move offsets after stim times
        stim_move_off_times = all_move_off_times(cell2mat(arrayfun(@(X) find(all_move_on_times>X,1,'first'),stim_move_on_times, 'UniformOutput', 0)));
        behaviour(animal_id).stim_move_off_times{day_index} = stim_move_off_times;
        
        % reaction times
        reaction_times = stim_move_on_times - stimOn_times';
        behaviour(animal_id).reaction_times{day_index} = reaction_times;
        
    end
    
    % find indices for days of reversal task
    protocol = 'AP_stimWheelLeftReverse';
    experiments = AP_find_experiments(animal,protocol);
    reversal_task_days = {experiments.day};
    reversal_task_days_mask = ismember(behaviour(animal_id).day, reversal_task_days);
    behaviour(animal_id).reversal_task_days_mask = reversal_task_days_mask;
    
    % find indices for days for original task
    original_task_days_mask = (~ismember(behaviour(animal_id).day, reversal_task_days))&(~ismember(behaviour(animal_id).day, muscimol_days));
    behaviour(animal_id).original_task_days_mask = original_task_days_mask;
 
    disp(['Done with ' animal])
end

disp('Done all')

% Save 
save('all_mice_behaviour.mat', 'behaviour', '-v7.3')
disp('Saved')


