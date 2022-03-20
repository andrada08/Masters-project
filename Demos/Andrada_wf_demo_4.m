%% Files to download

% I've copied some more files to
% \\zserver.cortexlab.net\Lab\Share\ajpeters\for_Andrada\widefield_alignment,
% copy that folder to your computer

addpath(genpath('C:\Users\Andrada\Documents\GitHub\npy-matlab'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\widefield'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\Lilrig'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\PupilDetection_DLC'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\AP_scripts_cortexlab'));


%% Standardizing data for combining

% Each widefield recording has a unique set of spatial components (U),
% so to combine data across recordings we need to put our data into a set
% of common spatial components. My common spatial components, which I call
% my "master U", were made from SVDing U's across a bunch of recordings to
% capture consistent spatial patterns across recordings. I think that math
% is dubious, but it doesn't really matter: your U just needs to capture
% the variance in the data, so it can be a little arbitrary. 

% Let's load in an example data set
animal = 'AP025';
day = '2017-09-28';
experiment = 1;
verbose = true;
AP_load_experiment;

% I added the master U into your share folder, so just point this variable
% to wherever you've stored it on your local computer: 
master_u_fn = 'D:\Andrada\Master project\widefield_alignment\U_master';
load(master_u_fn);

% After loading above, you'll have a variable in your workspace U_master

% Converting widefield data from an experiment into this master U space
% needs 2 steps:

% 1) align the data (the master U is aligned to the master retinotopy, so
% the current day's data needs to be aligned to the master retinotopy)
Udf_aligned = AP_align_widefield(Udf,animal,day);

% 2) we use the ChangeU function to give us a new set of temporal
% components (V) corresponding to our master U
fVdf_Umaster = ChangeU(Udf_aligned,fVdf,U_master);

% EXERCISE: convince yourself (hopefully) that reconstructed data with the
% experiment U's and the master U's are similar. Try reconstructing a bunch
% of frames using Udf/fVdf and U_master/fVdf_Umaster and scroll through
% them side-by-side, do they look similar?


example_frames = 500:600;

example_fluorescence_experiment = AP_svdFrameReconstruct(Udf,fVdf(:,example_frames));
example_fluorescence_master = AP_svdFrameReconstruct(Udf_aligned,fVdf_Umaster(:,example_frames));

AP_image_scroll(example_fluorescence_experiment);
axis image;

AP_image_scroll(example_fluorescence_experiment);
axis image;


% EXERCISE: try comparing the ROIs from the experiment/master U's to check
% how similar they are. In order to do this, we'll need to define an ROI,
% and use the same ROI for both reconstructions.
%
% First draw an ROI on top of this figure using the 'roipoly' function
avg_im_aligned = AP_align_widefield(avg_im,animal,day);
figure;
imagesc(avg_im_aligned);
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
axis image off

roi = roipoly();
%
% You should have an ROI mask now (mask is true inside ROI, false
% outside),you can now put this into AP_svd_roi - so far you've only used
% this function to draw ROIs, but it can also take pre-drawn ROIs
experiment_U_trace = AP_svd_roi(Udf_aligned,fVdf,[],[],roi);
master_U_trace = AP_svd_roi(U_master,fVdf_Umaster,[],[],roi);

% Plot those traces on top of each other - are they similar? If so, it
% means we can swap out the experiment U for the master U and it gives us
% functionally the same data.

figure; 
plot(experiment_U_trace, 'b')
hold on
plot(master_U_trace, 'r--')

% Let's quantify how similar they are: calculate the fraction of explained
% variance when moving from the master U to the experiment U.

% When we're combining data, we don't necessarily want to keep all possible
% components since it takes up more room and is slower to work with. How
% many components should we keep? Get ROI traces using different number of
% components (3rd dimension in the U's, 1st dimension in the V's) and
% calculating the explained variance. Make a plot with explained variance
% vs. number of components.

experiment_U_trace = AP_svd_roi(Udf_aligned,fVdf,[],[],roi);
explained_variance = nan(size(U_master, 3), 1);
for num_comp=1:size(U_master, 3)
    master_U_trace = AP_svd_roi(U_master(:,:,1:num_comp),fVdf_Umaster(1:num_comp,:,:),[],[],roi);
    diff_roi_traces = experiment_U_trace - master_U_trace;
    explained_variance(num_comp) = 1 - var(diff_roi_traces)/var(experiment_U_trace);
end

figure;
plot(1:size(U_master, 3), explained_variance, 'o')

% 10 components are enough??? so little!


% The raw data is dominated by large slow events, but we care about the
% faster deconvolved data. Make the explained variance vs. number of
% components plot above but this time using deconvolved data - is there a
% difference?


deconvolved_fVdf = AP_deconv_wf(fVdf);
experiment_U_trace = AP_svd_roi(Udf_aligned,deconvolved_fVdf,[],[],roi);

deconvolved_fVdf_Umaster = AP_deconv_wf(fVdf_Umaster);

deconvolved_explained_variance = nan(size(U_master, 3), 1);
for num_comp=1:size(U_master, 3)
    master_U_trace = AP_svd_roi(U_master(:,:,1:num_comp),deconvolved_fVdf_Umaster(1:num_comp,:,:),[],[],roi);
    diff_roi_traces = experiment_U_trace - master_U_trace;
    deconvolved_explained_variance(num_comp) = 1 - var(diff_roi_traces)/var(experiment_U_trace);

end

% takes too long!!!

figure;
plot(1:size(U_master, 3), deconvolved_explained_variance, 'o')

%% Combining data

% When we do our full analyses, we'll want to work on all of our data in a
% format that makes it easy to load and combine across recordings.

% EXERCISE: for a few of your mice (AP107/8/9 are probably a good bet),
% grab and save their passive data so you can analyze it all together. It
% will need these steps:
% 1) loop through each animal/recording
% 2) load that day's data
% 3) convert that day's data into the master U format
% 4) grab activity aligned to each stimulus
% 5) store the the activity and trial information, move to next loop
% 6) save

animals = {'AP107'} %,'AP108', 'AP109'};

master_u_fn = 'D:\Andrada\Master project\widefield_alignment\U_master';
load(master_u_fn);

master_activity = struct;

% create matrix of times for movie
timestep = 0.01;
start_time = -0.5;
end_time = 1;
timevec = [start_time:timestep:end_time];
stim_frame = (-start_time)*(1/timestep)+1;


master_activity.timestep = timestep;
master_activity.start_time = start_time;
master_activity.end_time = end_time;
master_activity.stim_frame = stim_frame;
master_activity.timevec = timevec;

for animal_id=1:length(animals)
    animal = animals{animal_id};
    master_activity(animal_id).animal = animal;
    protocol = 'AP_lcrGratingPassive';
    experiments = AP_find_experiments(animal,protocol);
    for day_index=1:length(experiments)
        day = experiments(day_index).day;
        master_activity(animal_id).day{day_index} = day;
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
        
        % temp ------------------------------------------------------------------------------------------------
%         
%         temp_all_stim_avg_act = nan(size(all_stim_act,1),size(all_stim_act,2),length(possible_stimuli));
%         % completed_trialStimulusValue = trialStimulusValue(1:n_trials);
%         
%         for stim_idx =1:length(possible_stimuli)
%             this_stim_act = all_stim_act(:,:,trialStimulusValue==possible_stimuli(stim_idx));
%             temp_all_stim_avg_act(:,:,stim_idx) = nanmean(this_stim_act,3);
%         end
%         
%         temp_all_stim_avg_act = temp_all_stim_avg_act - temp_all_stim_avg_act(:,stim_frame,:);
%         
%         deconvolved_all_stim_avg_act = AP_deconv_wf(temp_all_stim_avg_act, [], 1/timestep);
%  
%         % get fluoresence
%         all_stim_interval_avg_fluorescence = AP_svdFrameReconstruct(U_master,deconvolved_all_stim_avg_act);
%         
%         % video
%         AP_image_scroll(all_stim_interval_avg_fluorescence,timevec);
%         axis image;
        
        % save in struct
        master_activity(animal_id).stim_activity{day_index} = all_stim_act;
        master_activity(animal_id).trial_id{day_index} = trialStimulusValue;
    end
    disp(['Done with ' animal])
end

disp('Done all')

% Save 
save('master_activity.mat', 'master_activity', '-v7.3')
disp('Saved')


%% Load and plot
% Then write some code to load in the data and plot the average stimulus
% response across mice for each day (e.g. day 1 stimulus response averaged
% across all mice). What changes in the activity?

load('master_activity.mat');

master_u_fn = 'D:\Andrada\Master project\widefield_alignment\U_master';
load(master_u_fn);

animals = {master_activity.animal};

timestep = master_activity.timestep;
start_time = master_activity.start_time;
end_time = master_activity.end_time;
stim_frame = master_activity.stim_frame;
timevec = master_activity.timevec;

% do something to get the min/max number - though realistically mice all
% train on simple task for two weeks

num_days = length(master_activity(1).day);
possible_stimuli = [-1 0 1];

% to have dim for avg activity from size of first master activity
stim_act = master_activity(1).stim_activity{1};
all_stim_avg_act = nan(size(stim_act,1),size(stim_act,2),length(possible_stimuli));

for day_idx=1:num_days
    all_stim_avg_act = nan(size(stim_act,1),size(stim_act,2),length(possible_stimuli));
    for stim_idx =1:length(possible_stimuli)
        this_stim_act = [];
        for animal_id=1:length(animals)
            stim_act = master_activity(animal_id).stim_activity{day_idx};
            trialStimulusValue = master_activity(animal_id).trial_id{day_idx};
            this_stim_act = cat(3, this_stim_act, stim_act(:,:,trialStimulusValue==possible_stimuli(stim_idx)));
        end
        
        all_stim_avg_act(:,:,stim_idx) = nanmean(this_stim_act,3);    
    end
    all_stim_avg_act = all_stim_avg_act - all_stim_avg_act(:,stim_frame,:);
    
    deconvolved_all_stim_avg_act = AP_deconv_wf(all_stim_avg_act, [], 1/timestep);
    
    % get fluoresence
    all_stim_interval_avg_fluorescence = AP_svdFrameReconstruct(U_master,deconvolved_all_stim_avg_act);
    
    % video
    AP_image_scroll(all_stim_interval_avg_fluorescence,timevec);
    axis image;
end


% do videos for code where save for first day and compare to passive stim
% stuff - looks ok not the same but ok (because one is for master u and v)
% - with normal V and U 
% - then master ones 






