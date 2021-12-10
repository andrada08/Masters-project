%% Paths to stuff
addpath(genpath(cd));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\npy-matlab'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\widefield'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\AP_scripts_cortexlab'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\Lilrig'));

%% Move time difference

animal = 'AP109';
day = '2021-12-09';
experiment = 1;
verbose = true;
AP_load_experiment;

% use diff to get all movement onsets
t = Timeline.rawDAQTimestamps;

tmp_move = [0; diff(wheel_move)];
all_move_frames = find(tmp_move==1);
all_move_times = t(all_move_frames);

move_difference = diff(all_move_times);
mean_move_difference = mean(move_difference);

% prctile(move_difference,95)

% movement onset after stimulus
moveOn_times = all_move_times(cell2mat(arrayfun(@(X) find(all_move_times>X,1,'first'),stimOn_times', 'UniformOutput', 0)));
moveOn_frames = cell2mat(arrayfun(@(X) find(t==X,1),moveOn_times, 'UniformOutput', 0));

moveOn_difference = diff(moveOn_times);
mean_moveOn_difference = mean(moveOn_difference);

figure
subplot(2,2,1)
plot(moveOn_difference)
yline(mean_moveOn_difference)

subplot(2,2,2)
plot(move_difference)
yline(mean_move_difference)

title('AP109 2021-12-09')
hold on;

animal = 'AP109';
day = '2021-12-10';
experiment = 1;
verbose = true;
AP_load_experiment;

% use diff to get all movement onsets
t = Timeline.rawDAQTimestamps;

tmp_move = [0; diff(wheel_move)];
all_move_frames = find(tmp_move==1);
all_move_times = t(all_move_frames);

move_difference = diff(all_move_times);
mean_move_difference = mean(move_difference);

% prctile(move_difference,95)

% movement onset after stimulus
moveOn_times = all_move_times(cell2mat(arrayfun(@(X) find(all_move_times>X,1,'first'),stimOn_times', 'UniformOutput', 0)));
moveOn_frames = cell2mat(arrayfun(@(X) find(t==X,1),moveOn_times, 'UniformOutput', 0));

moveOn_difference = diff(moveOn_times);
mean_moveOn_difference = mean(moveOn_difference);

% plot
subplot(2,2,3)
plot(moveOn_difference)
yline(mean_moveOn_difference)

subplot(2,2,4)
plot(move_difference)
yline(mean_move_difference)

title('AP109 2021-12-10')


% AP108
animal = 'AP108';
day = '2021-12-09';
experiment = 1;
verbose = true;
AP_load_experiment;

% use diff to get all movement onsets
t = Timeline.rawDAQTimestamps;

tmp_move = [0; diff(wheel_move)];
all_move_frames = find(tmp_move==1);
all_move_times = t(all_move_frames);

move_difference = diff(all_move_times);
mean_move_difference = mean(move_difference);

% prctile(move_difference,95)

% movement onset after stimulus
moveOn_times = all_move_times(cell2mat(arrayfun(@(X) find(all_move_times>X,1,'first'),stimOn_times', 'UniformOutput', 0)));
moveOn_frames = cell2mat(arrayfun(@(X) find(t==X,1),moveOn_times, 'UniformOutput', 0));

moveOn_difference = diff(moveOn_times);
mean_moveOn_difference = mean(moveOn_difference);

% plot
figure
subplot(1,2,1)
plot(moveOn_difference)
yline(mean_moveOn_difference)

subplot(1,2,2)
plot(move_difference)
yline(mean_move_difference)

title('AP108')

%% Rightward move after reward for AP107

animal = 'AP107';
day = '2021-12-09';
experiment = 1;
verbose = true;
AP_load_experiment;

t = Timeline.rawDAQTimestamps;

tmp_move = [0; diff(wheel_move)];
all_move_frames = find(tmp_move==1);
all_move_times = t(all_move_frames);

move_reward_times = all_move_times(cell2mat(arrayfun(@(X) find(all_move_times>X,1,'first'),reward_t_timeline, 'UniformOutput', 0)));
move_reward_frames = cell2mat(arrayfun(@(X) find(t==X,1),move_reward_times, 'UniformOutput', 0));

figure
plot(wheel_velocity(move_reward_frames))


frames = 1:100;
after_move_frames = move_reward_frames + frames';
move_reward_velocity = mean(wheel_velocity(after_move_frames),1);

plot(move_reward_velocity)

% rightmove_frames = find(wheel_velocity(move_reward_frames)>0);



%% activity??

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
