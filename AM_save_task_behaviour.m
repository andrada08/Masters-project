%% Save task behaviour for AP107

animals = {'AP107'}; %,'AP108', 'AP109'};

behaviour = struct;

% create matrix of times for movie
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
    protocol = 'AP_stimWheel';
    flexible_name = true;
    experiments = AP_find_experiments(animal,protocol,flexible_name);
    for day_index=1:length(experiments)
        day = experiments(day_index).day;
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
        
    end
    disp(['Done with ' animal])
end

disp('Done all')

% Save 
save('AP107_behaviour.mat', 'behaviour', '-v7.3')
disp('Saved')


%% Load and plot reaction times 

% look through wf_demo_3 for this