%% Paths

addpath(genpath(cd));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\npy-matlab'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\widefield'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\AP_scripts_cortexlab'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\Lilrig'));

%% Load example dataset 
% (for the next things, let's use the same example dataset as before)

% Here's an example animal/day
animal = 'AP025';
day = '2017-09-28';
experiment = 1;
verbose = true;

AP_load_experiment;

%% Exploring data

% Here are two functions which could be useful for exploring data.

% Pixel correlation viewer: this plots the correlation of a particular
% pixel with all other pixels. This doesn't tell you much about an
% experiment, but I think it's useful for getting an intuition for brain
% areas. 
pixelCorrelationViewerSVD(Udf,fVdf);

% You can now click on any pixel, and it will show you the correlation
% between fluorescence in that pixel and all other pixels (yellow areas are
% highly correlated, black is 0 correlation, blue is negative
% correlation). Also - you can click on the figure and hit 'h' which puts
% it into 'hover mode', so now you can just move your mouse around to
% change which pixel is the reference pixel.

% EXERCISE: align the U's, deconvolve the V's, then feed those into the
% pixel correlation viewer. Then overlay the atlas regions on the image and
% look around at the data:

Udf_aligned = AP_align_widefield(Udf,animal,day);
fVdf_deconvolved = AP_deconv_wf(fVdf);

pixelCorrelationViewerSVD(Udf_aligned,fVdf_deconvolved);
AP_reference_outline('ccf_aligned','k');


% - Based on correlated areas, if you were to divide the brain into a few
% correlated chunks, where would the be? What areas/functions do those
% relate to?

% blobs each side - V1 and the secondary areas on that side
% also more lateral (temporal???) V1 has a correlation to the other side

% blobs on frontotemporal sides (primary motor cortices and top part of MOs and a somatosensory area)
% -- related to mouth/face movement

% blob frontomedial region (MOs, ALM etc.)

% blob dorsomedial? region (not sure what it is) - retrosplenial area (value, history, up/down tilt of head?)
% --- input from hippocampus

% mid-temporal blob on each side (somatosensory areas+motor)
% --- limb movement areas

% blob on barrel somatosensory area+bit of motor+MOs

% - Sometimes you'll see a pattern where you'll have separated spots of
% correlated pixels that merge together as you move the mouse. Where are
% these? Why do they show up where they do?

% moving from somatosensory towards frontal region - merges into that
% frontotemporal blob

% moving from middle (like bottom of MOs)/ from the the dorsomedial? region
% to a somatosensory area

% moving somatosensory + mid blob travelling up through MOs and MOp creates the frontotemporal blobs

% This next tool is a little more useful for exploring data as it relates
% to the experiment: it plots the face camera, eye camera, and widefield
% images (after you open, you can scroll with the scroll bar or the scroll
% wheel on the mouse).
AP_wfmovies(Udf,fVdf,frame_t,eyecam_fn,eyecam_t,facecam_fn,facecam_t)

% EXERCISE: look through the data and come up with some observations about
% the relationship between behavior and activity. For example: what areas
% are active while the mouse is moving? What's going on in the brain when
% the mouse isn't moving? Are there types of activity that you see often?
% Are there types of activity you see rarely? 

% when the mouse moves activity starts in the frontal regions and travels
% to the primary motor region and then moves to somatosensory??
% (and it travels towards lateral V1)
% before that there is activity in V1 on one side - the visual stimulus 

% after move somatosensory regions go dark

% licking?? - activity up moving to somatosensory

% mouse sitting still (960 onwards for example) - bit random??

% Then: run that tool using the deconvolved V's instead. You'll probably
% need to adjust the color balance as we did before - click on the
% widefield image, and in the command line type caxis(caxis/2) and repeat
% that command until you can see activity well. What's the difference
% between this and the non-deconvolved data? Can you make any observations
% about activity/behavior in this data that you couldn't see in the raw
% data?

AP_wfmovies(Udf,fVdf_deconvolved,frame_t,eyecam_fn,eyecam_t,facecam_fn,facecam_t)

% can see the areas better - activity doesn't linger so can see fluoresence
% travelling from a region to another

% activity in retrosplenial cortex

% oscillation in back when mouse is still - visual+retrosplenial cortex

%% Wheel movement

% The last component of our experiments is the wheel that the mice can
% turn, which is recorded into Timeline.

% EXERCISE: the wheel is recorded into Timeline under the name
% 'rotaryEncoder'. Plot the raw wheel trace in Timeline with time in
% seconds on the x-axis. This is the cumulative wheel position over the
% course of the whole experiment (i.e. each point is relative to the very
% start of the experiment, which is set at 0). Positive values means the
% wheel is turning right, negative values means the wheel is turning left.

% rotaryEncoder_index = strcmp({Timeline.hw.inputs.name}, 'rotaryEncoder');
% figure;
% plot(Timeline.rawDAQTimestamps,Timeline.rawDAQData(:,rotaryEncoder_index));
% xlabel('Time (s)');
% ylabel('Rotary Encoder');

% something is wrong here it doesn't look right - use wheel_position


figure;
plot(Timeline.rawDAQTimestamps,wheel_position);
xlabel('Time (s)');
ylabel('Rotary Encoder');

% EXERCISE: using the interp1 function, get and plot the wheel position
% aligned to each stimulus side/contrast. Note - you'll want "0" in your
% plot to be the wheel position before the stimulus comes on, so baseline
% subtraction will be important here to make sense of the data. Do these
% traces make sense given what the task is?

% find possible stimuli
possible_contrasts = unique(signals_events.trialContrastValues);
possible_stimuli = [-1;1].*possible_contrasts;
possible_stimuli = unique(possible_stimuli);

% find value of stimulus per trial
trialStimulusValue = signals_events.trialContrastValues .* signals_events.trialSideValues;

% create time matrix around stimulus onsets
timestep = 0.01;
start_time = -4;
end_time = 4;
timevec = [start_time:timestep:end_time];

stim_frame = (-start_time)*(1/timestep)+1;

time_stimulus = stimOn_times+timevec;

% find wheel position for stimulus onsets
all_stim_wheel_position = interp1(Timeline.rawDAQTimestamps,wheel_position,time_stimulus);
all_stim_wheel_position = permute(all_stim_wheel_position, [2,1]);

all_stim_avg_wheel_position = nan(size(all_stim_wheel_position,1),length(possible_stimuli));
completed_trialStimulusValue = trialStimulusValue(1:n_trials);

for stim_idx =1:length(possible_stimuli)
    this_stim_wheel = all_stim_wheel_position(:,completed_trialStimulusValue==possible_stimuli(stim_idx));
    all_stim_avg_wheel_position(:,stim_idx) = nanmean(this_stim_wheel,2);
end

% substract frame prior to stimulus
all_stim_avg_wheel_position = all_stim_avg_wheel_position - all_stim_avg_wheel_position(stim_frame-1,:);

% figure 
% set colourmap to differentiate left from right stimuli
map_left = [0 0 0.2
    0 0 0.4
    0 0 0.6
    0 0 0.8
    0 0 1.0];
map_zero = [0 1 0];
map_right = [0.2 0 0 
    0.4 0 0
    0.6 0 0 
    0.8 0 0 
    1.0 0 0];
map = [map_left; map_zero; map_right];

wheel_and_stim_figure = figure;
plot(all_stim_avg_wheel_position);
colororder(wheel_and_stim_figure,map);
xlabel('Frames')
xlim([1 length(timevec)])
ylabel('Wheel position')
legend(num2str(possible_stimuli));

% yes! when the stimulus is on the right the wheel is moved left and the
% opposite 


% EXERCISE: instead of using the wheel position, we can also look at wheel
% velocity. Plot the derivative of the wheel position using the 'diff'
% function, and zoom into the plot: you should be able to see discrete
% positive or negative pulses, which is what the wheel actually sends out
% to the computer.

calculated_wheel_velocity = [0; diff(wheel_position)];
figure;
plot(Timeline.rawDAQTimestamps, calculated_wheel_velocity);


% Those blips aren't super useful to look at, so I have a function to
% smooth these and break them into "movement" and "quiescence" periods.
% This is done in the loading script, so it's already available. Let's plot
% this processed velocity ('wheel_velocity'):
t = Timeline.rawDAQTimestamps;
figure; hold on;
plot(t,wheel_velocity,'k');
hold on;
% The movement/quiescence periods are stored in 'wheel_move': this is 0
% when the mouse is quiescent and 1 when the mouse is moving. Let's create
% a new trace which is the wheel velocity during movement and NaN during
% quiescence. NaNs (not-a-number) are placeholder values. In this case
% they're useful because NaNs aren't plotted, so we plot the wheel velocity
% during movement in a different color.
% (copy the wheel velocity trace)
wheel_velocity_movement = wheel_velocity;
% (whenever the mouse is quiescent, set those values to NaN)
wheel_velocity_movement(~wheel_move) = NaN;
% (plot the wheel velocity during movement in red)
plot(t,wheel_velocity_movement,'r');
% Zoom into this plot to see how the trace is either black+flat or
% red+wobbly.

% EXERCISE: using the 'wheel_move' trace, pull out times of movement
% onsets. Align the widefield data to these movement onsets and plot the
% average widefield fluorescence aligned to movement. What's going on in
% this movie?

% use diff to get all movement onsets
tmp_move = [0; diff(wheel_move)];
all_move_frames = find(tmp_move==1);
all_move_times = t(all_move_frames);

% get activity for window around all movement onsets
timestep = 0.1;
start_time = -2;
end_time = 2;
timevec = [start_time:timestep:end_time];

move_frame = (-start_time)*(1/timestep)+1;

time_move = all_move_times'+timevec;
 
all_move_act = interp1(frame_t,fVdf',time_move);
all_move_act = permute(all_move_act, [3,2,1]);

% get average across all moves and baseline substract 
all_move_avg_act = nanmean(all_move_act,3);
all_move_avg_act = all_move_avg_act - all_move_avg_act(:,move_frame,:);

all_move_avg_fluorescence = AP_svdFrameReconstruct(Udf,all_move_avg_act);

% movie 
AP_image_scroll(all_move_avg_fluorescence,timevec);
axis image;


% movement onset after stimulus
moveOn_times = all_move_times(cell2mat(arrayfun(@(X) find(all_move_times>X,1,'first'),stimOn_times', 'UniformOutput', 0)));
moveOn_frames = cell2mat(arrayfun(@(X) find(t==X,1),moveOn_times, 'UniformOutput', 0));

% get activity for window around all movement onsets
timestep = 0.1;
start_time = -2;
end_time = 2;
timevec = [start_time:timestep:end_time];

move_frame = (-start_time)*(1/timestep)+1;

time_move = moveOn_times'+timevec;
 
move_act = interp1(frame_t,fVdf',time_move);
move_act = permute(move_act, [3,2,1]);

% get average across all moves and baseline substract 
move_avg_act = nanmean(move_act,3);
move_avg_act = move_avg_act - move_avg_act(:,move_frame,:);

move_avg_fluorescence = AP_svdFrameReconstruct(Udf,move_avg_act);

% movie 
AP_image_scroll(move_avg_fluorescence,timevec);
axis image;

% all move onsets again
AP_image_scroll(all_move_avg_fluorescence,timevec);
axis image;

AP_image_scroll(all_move_avg_fluorescence-move_avg_fluorescence,timevec);
axis image;
colormap(brewermap([],'PRGn'));
caxis([-max(abs(caxis)),max(abs(caxis))])


% looks the same yay!!

% EXERCISE: is the movement-aligned activity different depending on how
% fast and in which direction the mouse is moving the wheel? This will be
% more complicated and will need a few steps: 
% - for each movement onset, get a number summarizing that movement (maybe
% the amplitude and sign of the velocity: remember positive values means
% the wheel is turning right, and negative values means it's turning left).
% One option for this is to split the velocity trace into a cell array,
% where each cell contains a full movement, which can be done with the
% mat2cell function. 
% - after you get a value for each movement representing how fast and in
% which direction the mouse is turning, you'll want to discretize those
% values into categories. Try using the 'discretize' function. 
% - once you have discrete values of movement (e.g. very fast left turn,
% small right turn, etc.), align the widefield images to those movement
% onsets and see how the activity varies with movement speed and direction.

% get mean of 10 frames? after
frames = 1:100;
after_move_frames = moveOn_frames + frames';
move_velocity = mean(wheel_velocity(after_move_frames),1);

% use prctile function - gives threshold for split - didn't really use this
move_percentiles = prctile(move_velocity,80);

% categorize values
move_categories = discretize(move_velocity,[-1 -0.04 -0.02 0 0.02 0.04 1],'categorical', ...
    {'quick left', 'medium left', 'slow left', 'slow right', 'medium right', 'quick right'}); 
possible_moves = unique(move_categories); 

move_avg_act = nan(size(move_act,1),size(move_act,2),length(possible_moves));

for move_idx =1:length(possible_moves)
    this_move_act = move_act(:,:,move_categories==possible_moves(move_idx));
    move_avg_act(:,:,move_idx) = nanmean(this_move_act,3);
end

move_avg_act = move_avg_act - move_avg_act(:,move_frame,:);

% get fluorescence 
move_interval_avg_fluorescence = AP_svdFrameReconstruct(Udf,move_avg_act);

% video
AP_image_scroll(move_interval_avg_fluorescence,timevec);
axis image;

figure;
plot(timevec,roi.trace')

% DEMO SEQUENTIAL COLORMAP
n = 5;
x = rand(100,n);
figure; hold on;
set(gca,'ColorOrder',copper(n));
plot(x)



%% Combining stimulus, movements, and imaging

% Now you've got the general tools to work with all of the data: let's try
% using them to get an example of separating stimulus- and movement-related
% activity. 

% EXERCISE: find the "reaction time" for each stimulus as the first time
% the mouse moves the wheel following the onset of each stimulus. As a
% sanity check, plot the median reaction time vs. stimulus and make sure
% this plot looks like what you expect it to look like. 

moveOn_times = moveOn_times';

all_reaction_times = nan(1,length(possible_stimuli));

for stim_idx =1:length(possible_stimuli)
    this_index = completed_trialStimulusValue==possible_stimuli(stim_idx);
    this_reaction_times = moveOn_times(this_index) - stimOn_times(this_index);
    all_reaction_times(stim_idx) = nanmean(this_reaction_times);
end

mean_reaction_time = nanmean(all_reaction_times);
median_reaction_time = median(all_reaction_times);

figure;
plot(possible_stimuli,all_reaction_times)
hold on;
yline(median_reaction_time)
yline(mean_reaction_time,'r')
ylabel('Reaction time (s)')
xlabel('Stimulus')

% for right there's much bigger values why

% EXERCISE: align widefield activity to the onset of all stimuli on the
% righthand screen (use all contrasts except zero - just lump them all
% together). Split these trials by reaction times into a few groups, so
% that you'll have trials with slow/medium/fast etc. reaction times. Make an
% average stimulus-aligned activity for each one of these reaction time
% groups. What changes in the activity across these reaction time groups?


rightstimOn_times = stimOn_times(completed_trialStimulusValue>0);

% get activity 
timestep = 0.01;
start_time = -2;
end_time = 2;
timevec = [start_time:timestep:end_time];

stim_frame = (-start_time)*(1/timestep)+1;

time_stimulus = rightstimOn_times+timevec;

% find activity for above time
right_stim_act = interp1(frame_t,fVdf',time_stimulus);
right_stim_act = permute(right_stim_act, [3,2,1]);

% get categories of reaction times
right_reaction_times = arrayfun(@(X) moveOn_times(find(moveOn_times>X,1,'first'))- X, rightstimOn_times, 'Uni', 1);

% right_reaction_times_percentiles = prctile(right_reaction_times,20)
% right_reaction_times_percentiles = prctile(right_reaction_times,40)
% right_reaction_times_percentiles = prctile(right_reaction_times,60)
% right_reaction_times_percentiles = prctile(right_reaction_times,80)
% right_reaction_times_percentiles = prctile(right_reaction_times,100)

right_reaction_times_categories = discretize(right_reaction_times,[0 1 2 4 6 22]); 
possible_categories = unique(right_reaction_times_categories);

right_stim_reaction_times_avg_act = nan(size(right_stim_act,1),size(right_stim_act,2),length(possible_categories));

for category=1:length(possible_categories)
    this_category_act = right_stim_act(:,:,right_reaction_times_categories==category);
    right_stim_reaction_times_avg_act(:,:,category) = nanmean(this_category_act,3);
end

right_stim_reaction_times_avg_act = right_stim_reaction_times_avg_act - right_stim_reaction_times_avg_act(:,stim_frame,:);

right_stim_reaction_times_avg_fluorescence = AP_svdFrameReconstruct(Udf,right_stim_reaction_times_avg_act);

% not sure
AP_image_scroll(right_stim_reaction_times_avg_fluorescence,timevec);
axis image;

% more frontal activity with fastest RT
% more lateralized visual response in contralateral side with slower reaction times




