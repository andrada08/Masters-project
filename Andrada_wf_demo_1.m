%% Repos to download

% Clone these repos: 
% https://github.com/petersaj/AP_scripts_cortexlab
% https://github.com/cortex-lab/widefield
% https://github.com/kwikteam/npy-matlab

% Do I do this??
addpath(genpath(cd));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\npy-matlab'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\widefield'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\AP_scripts_cortexlab'));

%% Load example dataset 

% Here's an example animal/day
animal = 'AP025';
day = '2017-09-28';
experiment = 1;
verbose = true;

% This is my catch-all code for loading data
% (finds the data on the server, loads, processes)
AP_load_experiment;

%% Timeline introduction
% Timeline is our code environment for input/output

% Anything in the experiment with a signal is recorded into timeline, which
% is saved into this structure:
Timeline

% The names of the inputs into timeline are stored here: 
{Timeline.hw.inputs.name};

% The timestamps for recorded signals are stored here: 
Timeline.rawDAQTimestamps;

% The recorded data is stored here (N timestamps x N signals)
Timeline.rawDAQData;

% The main inputs we'll use are: 
% photoDiode - this detects changes on the screen (e.g. stim presentation)
% rotaryEncoder - this is the position of the mouse's wheel
% rewardEcho - this is when a reward was given
% pcoExposure - this is when the widefield camera took a picture

% Example: if we want to plot the photodiode signal, we do that with this:
% 1) figure out which input belongs to the reward
photodiode_index = strcmp({Timeline.hw.inputs.name}, 'photoDiode');
% 2) use that to index which raw data to plot
figure;
plot(Timeline.rawDAQTimestamps,Timeline.rawDAQData(:,photodiode_index));
xlabel('Time (s)');
ylabel('Photodiode (volts)');

% EXERCISE: find the times that a reward was given by determining when the
% 'rewardEcho' signal turned on. That's been done in the load script, so
% compare your answer with 'reward_t_timeline'

reward_index = strcmp({Timeline.hw.inputs.name}, 'rewardEcho');
reward_data = Timeline.rawDAQData(:,reward_index);
reward_der = [0; diff(reward_data>reward_thresh)];
reward_times = Timeline.rawDAQTimestamps(find(reward_der==1));

plot(reward_times-reward_t_timeline)

% check 
all(reward_times == reward_t_timeline)

% plot
% figure;
% plot(Timeline.rawDAQTimestamps,reward_data);
% hold on 
% plot(reward_times,reward_data(reward),'*')
% xlabel('Time (s)');
% ylabel('Reward');

%% Signals introduction
% Signals is the code environment for our experiment protocols

% Signals saves experiment information for each trial, including things
% that can't be recorded by timeline (e.g. which stimulus was shown?)

% Each of these signals has both a value and a time, we usually only care
% about the values (because more accurate timings are in Timeline). These
% are stored in signals_events.xValues and signals_events.xTimes, for
% example trial numbers are stored in signals_events.trialNumValues and
% signals_events.trialNumTimes.

% Signals most important for this task are
% The side of the stim on each trial:
signals_events.trialSideValues; % (-1 is on the left, 1 is on the right)
% The contrast of the stim on each trial
signals_events.trialContrastValues;
% Correct (hit) = 1 or incorrect (miss) = 0 for each trial:
signals_events.hitValues;

% EXERCISE: using the above signals, make a plot showing the fraction of
% correct trials for each unique stimulus (side and contrast).

correct_index = find(signals_events.hitValues==1);
 
possible_contrasts = unique(signals_events.trialContrastValues);
possible_stimuli = [-1;1].*possible_contrasts;
possible_stimuli = unique(possible_stimuli);

trialStimulusValue = signals_events.trialContrastValues .* signals_events.trialSideValues;
[correct_fractions, stimuli] = grpstats(signals_events.hitValues, trialStimulusValue, {'mean', 'gname'});
these_stimuli = cellfun(@str2num, stimuli);

% still doesn't look right? why??
plot(these_stimuli, correct_fractions)
xlabel('Possible stimuli')
ylabel('Fractions')

% Plot psychometric curve 
% way it turned from hit/miss and left/right stimulus (do left or right turn) vs stimulus

turnValues = signals_events.trialSideValues;
turnValues(correct_index)= -signals_events.trialSideValues(correct_index);
left_turnValues = turnValues == -1;

[left_turn_fractions, stimuli] = grpstats(left_turnValues, trialStimulusValue, {'mean', 'gname'});
new_stimuli = cellfun(@str2num, stimuli);


% still doesn't look right? why??
plot(new_stimuli, left_turn_fractions)
xlabel('Possible stimuli')
ylabel('Fractions')
 
% use accumarray once + unique functions second output to do psychometric curve

%% Widefield data introduction
% Widefield data is saved in SVD format rather than in pixels because it
% takes up way less space but contains most of the information. 

% This has two parts to it:
Udf; % U's: the spatial components (Y pixels x X pixels x N components)
fVdf; % V's: the temporal components (N components x N timepoints)
frame_t; % these are the timestamps for the temporal components (seconds)

% The spatial components represent modes that all the pixels can vary in.
% They decrease in how much variance they explain, so the first one looks
% like a brain and explains most of the variance, and the last one looks
% like noise and explains very little: 
figure; 
subplot(1,2,1);
imagesc(Udf(:,:,1));
axis image off;
title('Spatial component 1')

subplot(1,2,2);
imagesc(Udf(:,:,2000));
axis image off;
title('Spatial component 2000')

% The temporal components represent how strong each spatial component is at
% each timepoint. The first one is usually big (that component explains
% lots of variance and the last one is usually very small (that component
% explains very little variance) - try zooming in to see the variance in
% the last component:
figure; hold on
plot(frame_t,fVdf([1,end],:)');
xlabel('Time (s)');
ylabel('Component weight');
legend({'Temporal component 1','Temporal component 2'});

% Pixel values for each frame ('reconstructing') can then be gotten by 
% matrix multiplcation of the spatial and temporal components 
% (this has some reshapes: the spatial components need to be changed from 
% 3D to 2D for matrix multiplcation then the resulting frame has to be 
% changed from one long vector into a Y pixels x X pixels size)
example_frame = 500;
example_fluorescence_long = ...
    reshape(Udf,[],size(Udf,3))*fVdf(:,example_frame);
example_fluorescence = ...
    reshape(example_fluorescence_long,size(Udf,1),size(Udf,2));
figure;
imagesc(example_fluorescence);
axis image off;
title('Example frame fluorescence');

% We've got this function to do the above quickly, so for example we can
% reconstruct a few frames together, which makes a 3D matrix of Y pixels x
% X pixels x N frames
example_frames = 500:600;
example_fluorescence = AP_svdFrameReconstruct(Udf,fVdf(:,example_frames));
% and then we can view that with a function I have for scrolling through 3D
% matricies: 
AP_image_scroll(example_fluorescence);
axis image;

% The other nice thing about working with SVD data is that any linear
% operation can be done on the temporal component ('V space') before
% reconstruction ('pixel space'). For example, we can get the average
% activity across the whole day simply by taking the average V, then
% reconstructing: 
avg_V = nanmean(fVdf,2);
avg_fluorescence = AP_svdFrameReconstruct(Udf,avg_V);
figure;imagesc(avg_fluorescence);
axis image off
title('Average fluoresence');

% Most of the things we'd want to do are linear
% (adding/subtracting/multiplying/dividing: things that can be done in any
% order). Some things we might be interested in aren't linear so we'd have
% to do in pixel space (e.g. maximum fluorescence in each pixel, standard
% devation - because that involves a square root).

% EXERCISE: make a reward-triggered average movie, showing the average
% fluorescence -1:1 second around rewards. You found the reward times
% above, and since this is averaging you can do it in V-space before
% reconstructing into pixel space.

% Baseline-subtract the triggered average movie

% get all start and end indexes for the window
start_indices = [];
reward_indices = [];
end_indices = [];
for i = reward_times
    start_indices = [start_indices find(frame_t>=i-1,1,'first')];
    reward_indices = [reward_indices find(frame_t>=i,1,'first')];
    end_indices = [end_indices find(frame_t>=i+1,1,'first')];
end

% how many frames here
nframes = length(start_indices(1):end_indices(1));
reward_baseline_index = reward_indices(1)-start_indices(1)+1;

% get the matrix of activity and average across how many indices
reward_interval_act = cell2mat(arrayfun (@(X) fVdf(:,start_indices(X):end_indices(X)), [1:length(start_indices)], 'UniformOutput', 0));
reward_interval_act = reshape(reward_interval_act(:), size(fVdf,1), nframes, length(start_indices));
reward_interval_avg_V = nanmean(reward_interval_act,3);
reward_baseline = reward_interval_avg_V(:,reward_baseline_index);
reward_interval_avg_V = reward_interval_avg_V - reward_baseline; 
reward_interval_avg_fluorescence = AP_svdFrameReconstruct(Udf,reward_interval_avg_V);

AP_image_scroll(reward_interval_avg_fluorescence);
axis image;


% Make a triggered average for all stimuli

% get all start and end indexes for the window
start_indices = [];
stimulus_indices = [];
end_indices = [];
for i = stimOn_times'
    start_indices = [start_indices find(frame_t>=i-1,1,'first')];
    stimulus_indices = [stimulus_indices find(frame_t>=i,1,'first')];
    end_indices = [end_indices find(frame_t>=i+1,1,'first')];
end

% how many frames here
nframes = length(start_indices(1):end_indices(1));
stimulus_baseline_index = stimulus_indices(1)-start_indices(1)+1;

% get the matrix of activity and average across how many indices
stimulus_interval_act = cell2mat(arrayfun (@(X) fVdf(:,start_indices(X):end_indices(X)), [1:length(start_indices)], 'UniformOutput', 0));
stimulus_interval_act = reshape(stimulus_interval_act(:), size(fVdf,1), nframes, length(start_indices));
stimulus_interval_avg_V = nanmean(stimulus_interval_act,3);
stimulus_baseline = stimulus_interval_avg_V(:,stimulus_baseline_index);
stimulus_interval_avg_V = stimulus_interval_avg_V - stimulus_baseline; 
stimulus_interval_avg_fluorescence = AP_svdFrameReconstruct(Udf,stimulus_interval_avg_V);

AP_image_scroll(stimulus_interval_avg_fluorescence);
axis image;


% split by stimuli 
% interp1 function and sampling rate is 35Hz but change timesteps to see
% what happens - timevec with timesteps in it for 1 second around stimulus
% onset and then use function and then can put a 3D matrix into the scroll
% thing!!!!

timestep = 0.1;
timevec = [-1:timestep:1];

time_stimulus = stimOn_times+timevec;
 
all_stimuli_act = interp1(frame_t,fVdf',time_stimulus);

avg_check = permute(nanmean(all_stimuli_act,1), [3, 2, 1]);
fluorescence_check = AP_svdFrameReconstruct(Udf,avg_check);

AP_image_scroll(fluorescence_check,timevec);
axis image;

% split by stimuli use only one for loop!!!

trialStimulusValue 






