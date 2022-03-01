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

explained_variance = nan(size(U_master, 3), 1);
for num_comp=1:size(U_master, 3)
    master_U_trace = AP_svd_roi(U_master(:,:,1:num_comp),deconvolved_fVdf_Umaster(1:num_comp,:,:),[],[],roi);
    diff_roi_traces = experiment_U_trace - master_U_trace;
    explained_variance(num_comp) = 1 - var(diff_roi_traces)/var(experiment_U_trace);

end

% takes too long!!!

figure;
plot(1:size(U_master, 3), explained_variance, 'o')

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
%
% Then write some code to load in the data and plot the average stimulus
% response across mice for each day (e.g. day 1 stimulus response averaged
% across all mice). What changes in the activity?

animals = {'AP107','AP108'};
% 
% for animal_id=1:length(animals)
%     animal = animals{animal_id};
%     eyeblink(animal_id).animal = animal;
%     protocol = 'AP_lcrGratingPassive';
%     experiments = AP_find_experiments(animal,protocol);
%     eyeblink(animal_id).blinks = nan(length(experiments),3);
%     for day_index=1:length(experiments)
%         day = experiments(day_index).day;
%         eyeblink(animal_id).day{day_index} = day;
%         experiment = experiments(day_index).experiment(end);
%         load_parts.imaging = false;
%         load_parts.cam = true;
%         verbose = true;
%         AP_load_experiment;
%         












