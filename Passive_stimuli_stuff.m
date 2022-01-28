%% Paths to stuff
addpath(genpath(cd));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\npy-matlab'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\widefield'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\AP_scripts_cortexlab'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\Lilrig'));

%% Initial stuff
% mice = {'AP107','AP108','AP109'};
% all_dates = {'2021-11-23','2021-11-24','2021-12-07','2021-12-08','2021-12-09','2021-12-10','2021-12-11','2021-12-12','2021-12-13','2021-12-14','2021-12-15'};
% passive_only_dates = {'2021-11-23','2021-11-24','2021-12-07'};
% training_dates = {'2021-12-08','2021-12-09','2021-12-10','2021-12-11','2021-12-12','2021-12-13','2021-12-14','2021-12-15'};
% dates = {'2021-12-07'};

mice = {'AP110', 'AP111', 'AP112'};
passive_only_dates = {'2022-01-19','2022-01-20','2022-01-21'};

% retinotopy for each mouse
% mouse_idx = 1;
% 
% for mouse_idx=1:length(mice)
%     animal = mice{mouse_idx};
%     day = '2021-12-07';
%     experiment = 2;
%     verbose = true;
%     AP_load_experiment;
%     
%     lilrig_retinotopy
% end

%% Load experiment
% for day_idx=1:length(training_dates)
%     for mouse_idx=1:length(mice)
%         animal = mice{mouse_idx};
%         day = training_dates{day_idx};
%         experiment = 2;
%         if strcmp(day,'2021-12-13') && strcmp(animal,'AP107')
%             experiment = 3; 
%         end
% %         if day_idx<=3
% %             experiment = 1;
% %         end

for day_idx=1:length(passive_only_dates)
    for mouse_idx=1:length(mice)
        animal = mice{mouse_idx};
        day = passive_only_dates{day_idx};
        experiment = 1;
        verbose = true;
        AP_load_experiment;
        
        %% Processing stuff
        
        % alignment and deconvolution???
        % test_aligned = AP_align_widefield(Udf,animal,day);
        % deconvolved_fVdf = AP_deconv_wf(fVdf);
        
        %% Stimulus movies
        
        
        % find stimuli
        % trialStimulusValue = trial_conditions(:,2)./90;
        % possible_stimuli = unique(trialStimulusValue);
        
        
        % find possible stimuli
        possible_contrasts = unique(signals_events.stimContrastValues);
        possible_sides = unique(signals_events.stimAzimuthValues)/90;
        possible_stimuli = possible_sides.*possible_contrasts;
        possible_stimuli = unique(possible_stimuli);
        
        % find value of stimulus per trial
        trialStimulusValue = signals_events.stimAzimuthValues/90 .* signals_events.stimContrastValues;
        
        % create matrix of times for movie
        timestep = 0.01;
        start_time = -2;
        end_time = 3;
        timevec = [start_time:timestep:end_time];
        
        stim_frame = (-start_time)*(1/timestep)+1;
        
        time_stimulus = stimOn_times+timevec;
        
        % find activity for above time
        all_stim_act = interp1(frame_t,fVdf',time_stimulus);
        all_stim_act = permute(all_stim_act, [3,2,1]);
        
        % calculate average act for each stimulus
        all_stim_avg_act = nan(size(all_stim_act,1),size(all_stim_act,2),length(possible_stimuli));
        % completed_trialStimulusValue = trialStimulusValue(1:n_trials);
        
        for stim_idx =1:length(possible_stimuli)
            this_stim_act = all_stim_act(:,:,trialStimulusValue==possible_stimuli(stim_idx));
            all_stim_avg_act(:,:,stim_idx) = nanmean(this_stim_act,3);
        end
        
        all_stim_avg_act = all_stim_avg_act - all_stim_avg_act(:,stim_frame,:);
        
        % get fluoresence
        all_stim_interval_avg_fluorescence = AP_svdFrameReconstruct(Udf,all_stim_avg_act);
        
        % video
        AP_image_scroll(all_stim_interval_avg_fluorescence,timevec);
        axis image;
    end
end