%% Repos to download 

% Clone these repos and add them to your path:
% https://github.com/cortex-lab/Lilrig

% I've copied some files to
% \\zserver.cortexlab.net\Lab\Share\ajpeters\for_Andrada\widefield_alignment,
% copy that folder to your computer

% You'll have to edit 2 functions to use the path where you copied the
% above files: 
% edit AP_align_widefield, enter the path in line 30
% edit AP_reference_outline, enter the path in line 30

addpath(genpath(cd));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\npy-matlab'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\widefield'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\AP_scripts_cortexlab'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\Lilrig'));

%% Retinotopy, alignment, and atlas overlay
% Now that we can work with widefield data, we want to know where different
% brain areas are in our image. The atlas we use is from the Allen
% Institute, sometimes called the common coordinate framework (CCF) or the
% Allen reference atlas (ARA). Here's the paper describing it that has some
% nice figures listing all the top-down areas we see in widefield:
% https://www.sciencedirect.com/science/article/pii/S0092867420304025?via%3Dihub
%
% We need to align the CCF atlas to our image to figure out which areas are
% which, and we also need something to align from one animal to another so
% we can combine our data. We do this both of these with retinotopy.
%
% 'Retinotopy' means there's a 1:1 map of visual space onto the brain,
% which happens in all visual cortical areas - this gives us a few maps,
% one for primary visual cortex and one for all secondary visual areas. The
% most recent mapping we use for reference is in this paper: 
% https://elifesciences.org/articles/18372
%
% Each visual area in that paper has a color, red or blue, which is the
% direction of the map. In V1 the map is flipped from the screen: 
% leftward on the screen gives rightward on the brain, which is colored 
% as blue. This reverses at each visual step, so in visual area AM for 
% example leftward on the screen is leftward on the brain, which is 
% colored red. 
%
% We get our retinotopic maps from a "sparse noise" stimulus: show a bunch
% of dots on the screen at random, get a screen-to-brain mapping, then
% color each map by whether it is oriented the same way as the screen or
% flipped. Here's an example recording and retinotopic map:

% Load in a sparse noise experiment
animal = 'AP048';
day = '2019-09-30';
experiment = 1;
verbose = true;
AP_load_experiment;

% Run code to pull out the retinotopy (takes a few minutes):
lilrig_retinotopy

% That gives you the retinotopic map (blue/red) by itself, overlaid on the
% average image, and with the CCF atlas areas aligned and overlaid.

% I align mice by aligning the retinotopic map for each mouse to a "master"
% retinotopic map (which is the aligned average of a bunch of mice). I
% then store the alignment transformation for each mouse/day, which I load
% to do the alignment. In the first part of this script, you copied the
% alignment file and changed the paths in the relevant files. 

% Let's align the example day's retinotopic map, which is in the variable
% 'vfs_boot_mean': 
vfs_boot_mean_aligned = AP_align_widefield(vfs_boot_mean,animal,day);

% We can plot the aligned retinotopic map, and then use the function
% AP_reference_outline to draw the aligned CCF regions on top
figure; imagesc(vfs_boot_mean_aligned); axis image off;
AP_reference_outline('ccf_aligned','k');

% The big yellow blob is V1, which should fit into the top part of V1 in
% the CCF. The little yellow blob sticking out of the top of V1 is visual
% area AM - that should be in the corner of the triangle that's VisAM.

% For reference, let's do that plot for the unaligned data: you can see it
% doesn't match up as well (it's not super far off, because each animal is
% imaged in roughly the same orientation).
figure; imagesc(vfs_boot_mean); axis image off;
AP_reference_outline('ccf_aligned','k');

%% Load example dataset 
% (for the next things, let's use the same example dataset as before)

% Here's an example animal/day
animal = 'AP025';
day = '2017-09-28';
experiment = 1;
verbose = true;

% This is my catch-all code for loading data
% (finds the data on the server, loads, processes)
AP_load_experiment;

%% Align stimulus responses

% EXERCISE: you previously got the average widefield response to each
% stimulus - do that again, but with master-aligned data. Our widefield
% data comes in U's (spatial components) and V's (spatial components), so
% all you have to do is align the U's once (feed them into
% AP_align_widefield) to get aligned U's and then do your analysis on that
% data. View those movies with AP_image_scroll, and overlay the CCF areas
% with AP_reference_outline. 

% align spatial components ????
test_aligned = AP_align_widefield(Udf,animal,day);
% 
% AP_image_scroll(Udf);
% axis image;
% AP_image_scroll(test_aligned);
% axis image;

% imagesc(test_aligned)

% % test stuff
% %
% % imagesc(new_vfs_aligned)
% % 
% % figure
% % imagesc(vfs_boot_mean)
% 
% % make new Udf???
% imagesc(Udf(:,:,2000))

% find possible stimuli
possible_contrasts = unique(signals_events.trialContrastValues);
possible_stimuli = [-1;1].*possible_contrasts;
possible_stimuli = unique(possible_stimuli);

% find value of stimulus per trial
trialStimulusValue = signals_events.trialContrastValues .* signals_events.trialSideValues;

% create matrix of times for movie
timestep = 0.01;
start_time = -4;
end_time = 4;
timevec = [start_time:timestep:end_time];

stim_frame = (-start_time)*(1/timestep)+1;

time_stimulus = stimOn_times+timevec;

% find activity for above time
all_stim_act = interp1(frame_t,fVdf',time_stimulus);
all_stim_act = permute(all_stim_act, [3,2,1]);

% avg_check = permute(nanmean(all_stimuli_act,1), [3, 2, 1]);
% fluorescence_check = AP_svdFrameReconstruct(Udf,avg_check);
% 
% AP_image_scroll(fluorescence_check,timevec);
% axis image;

% split by stimuli use only one for loop!!!

all_stim_avg_act = nan(size(all_stim_act,1),size(all_stim_act,2),length(possible_stimuli));
completed_trialStimulusValue = trialStimulusValue(1:n_trials);

for stim_idx =1:length(possible_stimuli)
    this_stim_act = all_stim_act(:,:,completed_trialStimulusValue==possible_stimuli(stim_idx));
    all_stim_avg_act(:,:,stim_idx) = nanmean(this_stim_act,3);
end

all_stim_avg_act = all_stim_avg_act - all_stim_avg_act(:,stim_frame,:);

% use aligned U's instead of Udf
all_stim_interval_avg_fluorescence = AP_svdFrameReconstruct(test_aligned,all_stim_avg_act);

% all_stim_interval_avg_fluorescence = AP_svdFrameReconstruct(Udf,all_stim_avg_act);


% video
AP_image_scroll(all_stim_interval_avg_fluorescence,timevec);
axis image;
AP_reference_outline('ccf_aligned','k');

%% Regions of interest (ROIs)

% When we want to plot a line trace for fluorescence in one area, we can
% draw a region of interest (ROI), average all the pixels in the ROI, and
% plot the ROI fluorescence across time. 

% One way to do that is with the AP_image_scroll GUI: if you hit 'r', you
% can draw an ROI, then double click when you're done and it'll average
% pixels inside the ROI and save as the structure 'roi'. One field of 'roi'
% is 'roi.trace': that's N conditions x N timepoints.

% EXERCISE: plot the timecourse of stimulus-triggered activity in the
% primary visual cortex from the movies you made earlier. If you're looking
% in the left visual cortex, you should see more activity with increasing
% contrast on the right-hand screen.

AP_image_scroll(all_stim_interval_avg_fluorescence,timevec);
axis image;
AP_reference_outline('ccf_aligned','k');

saved_roi = roi;

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

% plot - red is right, blue is left
roi_figure = figure;
plot(saved_roi.trace');
colororder(roi_figure,map);
xlim([1,length(timevec)]);
xline(stim_frame);
legend(num2str(possible_stimuli));


% Another way to get activity in an ROI is with SVD data, with the function
% AP_svd_roi, the inputs are the U's, the V's, and a guide image. This
% returns the activity within an ROI across the whole recording. 

% EXERCISE: use the function below to get activity within an ROI over the
% primary visual cortex. Plot that over time along with markers (like
% lines) when the stimulus was presented: do you see a fluorescence
% response to the stimulus? 
roi_trace = AP_svd_roi(Udf,fVdf,avg_im);

figure
plot(frame_t,roi_trace)
xline(stimOn_times(completed_trialStimulusValue==1),'r')
% if zoom in see response

% EXERCISE: use interp1 to align and average the ROI to get average
% fluorescence for each stimulus presentation. This should look the same as
% the traces above where you made an ROI directly on the average movie -
% does it?

roi_stim_act = interp1(frame_t,roi_trace,time_stimulus);
roi_stim_act = permute(roi_stim_act, [2 1]);

roi_stim_avg_act = nan(size(roi_stim_act,1),length(possible_stimuli));
completed_trialStimulusValue = trialStimulusValue(1:n_trials);

for stim_idx =1:length(possible_stimuli)
    this_stim_act = roi_stim_act(:,completed_trialStimulusValue==possible_stimuli(stim_idx));
    roi_stim_avg_act(:,stim_idx) = nanmean(this_stim_act,2);
end

roi_stim_avg_act = roi_stim_avg_act - roi_stim_avg_act(stim_frame,:);

new_roi_figure = figure;
plot(roi_stim_avg_act);
colororder(new_roi_figure,map);
xlim([1,length(timevec)]);
xline(stim_frame);
legend(num2str(possible_stimuli));
title('ROI on activity then average');
ylim([-0.01 0.06]);

% other figure to compare
roi_figure = figure;
plot(saved_roi.trace');
colororder(roi_figure,map);
xlim([1,length(timevec)]);
xline(stim_frame);
legend(num2str(possible_stimuli));
title('ROI on average activity');
ylim([-0.01 0.06]);

%% Deconvolution 

% The raw fluorescence is much slower than spikes, and also has an
% expotential decay when there are no spikes. We can try to get rid of this
% by "deconvolving" our fluorescence - this should be reversing the GCaMP
% signal being the spike train convolved with an expotential decay. 

% The deconvolution kernel comes from simultaneous electrophysiology and
% imaging in the cortex, and tries to find a best match between the two.

% The function is 'AP_deconv_wf' and is used like 
% deconvolved_activity = AP_deconv_wf(raw_activity) where 'raw_activity' is
% an N-dimensional matrix there the second dimension is always time (e.g.
% could be conditions x time, or ROIs x time x conditions)

% EXERCISE: use AP_deconv_wf to plot a deconvolved stimulus-aligned
% fluorescence in V1. What's the difference between deconvolved and raw?

deconvolved_fVdf = AP_deconv_wf(fVdf);

deconvolved_roi_trace = AP_svd_roi(Udf,deconvolved_fVdf,avg_im);
plot(deconvolved_roi_trace)

deconvolved_roi_stim_act = interp1(frame_t,deconvolved_roi_trace,time_stimulus);
deconvolved_roi_stim_act = permute(deconvolved_roi_stim_act, [2 1]);

deconvolved_roi_stim_avg_act = nan(size(deconvolved_roi_stim_act,1),length(possible_stimuli));

for stim_idx =1:length(possible_stimuli)
    this_stim_act = deconvolved_roi_stim_act(:,completed_trialStimulusValue==possible_stimuli(stim_idx));
    deconvolved_roi_stim_avg_act(:,stim_idx) = nanmean(this_stim_act,2);
end

deconvolved_roi_stim_avg_act = deconvolved_roi_stim_avg_act - deconvolved_roi_stim_avg_act(stim_frame,:);

deconvolved_roi_figure = figure;
plot(deconvolved_roi_stim_avg_act);
colororder(deconvolved_roi_figure,map);
xline(stim_frame);
xlim([1,length(timevec)]);
legend(num2str(possible_stimuli));
title('Deconvolved data');

% raw version

raw_roi_figure = figure;
plot(roi_stim_avg_act);
colororder(raw_roi_figure,map);
xline(stim_frame);
xlim([1,length(timevec)]);
legend(num2str(possible_stimuli));
title('Raw data');


% EXERCISE: you can use AP_deconv_wf on the entire day with
% AP_deconv_wf(fVdf) (we only need to deconvolve in time - V's, we don't
% convolve in space so the U's stay the same). As before, make movies of
% average stimulus-aligned fluorescence, but using deconvolved data. What's
% the difference between deconvolved and raw?

deconvolved_fVdf = AP_deconv_wf(fVdf);

deconvolved_all_stim_act = interp1(frame_t,deconvolved_fVdf',time_stimulus);
deconvolved_all_stim_act = permute(deconvolved_all_stim_act, [3,2,1]);

deconvolved_all_stim_avg_act = nan(size(deconvolved_all_stim_act,1),size(deconvolved_all_stim_act,2),length(possible_stimuli));

for stim_idx =1:length(possible_stimuli)
    this_stim_act = deconvolved_all_stim_act(:,:,completed_trialStimulusValue==possible_stimuli(stim_idx));
    deconvolved_all_stim_avg_act(:,:,stim_idx) = nanmean(this_stim_act,3);
end

deconvolved_all_stim_avg_act = deconvolved_all_stim_avg_act - deconvolved_all_stim_avg_act(:,stim_frame,:);

deconvolved_all_stim_interval_avg_fluorescence = AP_svdFrameReconstruct(Udf,deconvolved_all_stim_avg_act);

% deconvolved
AP_image_scroll(deconvolved_all_stim_interval_avg_fluorescence,timevec);
axis image;

% raw
AP_image_scroll(all_stim_interval_avg_fluorescence,timevec);
axis image;

% activity doesn't persist in deconvolved as in raw


%% Hemodynamics

% Our imaging is done with fluorescence, excitation light of one wavelength
% (blue or violet) triggers light emission from GCaMP in a different
% wavelength (green). We use two different excitation colors because they
% each excite GCaMP in a different way: blue excites GCaMP differently
% depending on whether it's bound to calcium (bright with calcium, dark
% without calcium), and violet excites GCaMP the same regardless of whether
% it's bound to calcium. 

% The reason this is useful is to counteract hemodynamic effects: when an
% area is active, blood comes into that area, which occludes our GCaMP
% fluorescence. We want to have a separate measure of these hemodynamic
% effects (how much light is the blood blocking?) so that we can 'correct'
% our activity for any effects from blood that we don't want to measure.

% This is corrected in the loading script, but we can still look at these
% components separately. So far, you've been using the
% hemodynamically-corrected data (Udf, fVdf), where the uncorrected
% versions are in:
% blue (activity): Un, fVn (n for 'neural')
% violet (hemodynamic): Uh, fVh (h for 'hemodynamic');

% EXERCISE: Make a stimulus-triggered average movie of fluorescence from
% blue light (activity) and from violet light (hemodynamics). Plot
% stimulus-triggered average activity in V1 for both activity and
% hemodynamics. What are the main features of these plots, and why do they
% happen?


% blue - activity

activity_roi_trace = AP_svd_roi(Un,fVn,avg_im);

activity_roi_stim_act = interp1(frame_t,activity_roi_trace,time_stimulus);
activity_roi_stim_act = permute(activity_roi_stim_act, [2 1]);

activity_roi_stim_avg_act = nan(size(activity_roi_stim_act,1),length(possible_stimuli));

for stim_idx =1:length(possible_stimuli)
    this_stim_act = activity_roi_stim_act(:,completed_trialStimulusValue==possible_stimuli(stim_idx));
    activity_roi_stim_avg_act(:,stim_idx) = nanmean(this_stim_act,2);
end

activity_roi_stim_avg_act = activity_roi_stim_avg_act - activity_roi_stim_avg_act(stim_frame,:);

activity_roi_figure = figure;
plot(activity_roi_stim_avg_act);
colororder(activity_roi_figure,map);
xlim([1,length(timevec)]);
xline(stim_frame);
legend(num2str(possible_stimuli));
title('Activity');

% red - hemodynamics

hemodynamics_roi_trace = AP_svd_roi(Uh,fVh,avg_im);

hemodynamics_roi_stim_act = interp1(frame_t,hemodynamics_roi_trace,time_stimulus);
hemodynamics_roi_stim_act = permute(hemodynamics_roi_stim_act, [2 1]);

hemodynamics_roi_stim_avg_act = nan(size(hemodynamics_roi_stim_act,1),length(possible_stimuli));

for stim_idx =1:length(possible_stimuli)
    this_stim_act = hemodynamics_roi_stim_act(:,completed_trialStimulusValue==possible_stimuli(stim_idx));
    hemodynamics_roi_stim_avg_act(:,stim_idx) = nanmean(this_stim_act,2);
end

hemodynamics_roi_stim_avg_act = hemodynamics_roi_stim_avg_act - hemodynamics_roi_stim_avg_act(stim_frame,:);

hemodynamics_roi_figure = figure;
plot(hemodynamics_roi_stim_avg_act);
colororder(hemodynamics_roi_figure,map);
xlim([1,length(timevec)]);
xline(stim_frame);
legend(num2str(possible_stimuli));
title('Hemodynamics');

















