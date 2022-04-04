%% Paths to stuff
addpath(genpath(cd));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\npy-matlab'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\widefield'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\Lilrig'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\PupilDetection_DLC'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\AP_scripts_cortexlab'));

% animals = {'AP100','AP101','AP103','AP104','AP105','AP106','AP107','AP108','AP109'};

animals = {'AP107','AP108', 'AP109', 'AP110', 'AP111', 'AP112', 'AP113', 'AP114', 'AP115'};

eyeblink = struct;

for animal_id=1:length(animals)
    animal = animals{animal_id};
    eyeblink(animal_id).animal = animal;
    protocol = 'AP_lcrGratingPassive';
    experiments = AP_find_experiments(animal,protocol);
    eyeblink(animal_id).blinks = nan(length(experiments),3);
    for day_index=1:length(experiments)
        day = experiments(day_index).day;
        eyeblink(animal_id).day{day_index} = day;
        experiment = experiments(day_index).experiment(end);
        load_parts.imaging = false;
        load_parts.cam = true;
        verbose = true;
        AP_load_experiment;
        
        %% Load DLC output old
        
        if isempty(eyecam_dlc_filename)
            disp(['Eyecam DLC missing for ' animal ' ' day ' ' num2str(experiment)])
            continue
        end
        
        eyecam_DLC_output = readmatrix(eyecam_dlc_filename);
        
        pupil_top = eyecam_DLC_output(:,2:4); % [x y likelihood]
        pupil_bottom = eyecam_DLC_output(:,5:7);
        pupil_right = eyecam_DLC_output(:,8:10);
        pupil_left = eyecam_DLC_output(:,11:13);
        lid_top = eyecam_DLC_output(:,14:16);
        lid_bottom = eyecam_DLC_output(:,17:19);
        data = [pupil_top, pupil_bottom, pupil_left, pupil_right, lid_top, lid_bottom];
        
        %% Parameters - from pupil DLC page
        eyecam_DLC_params.minCertainty = 0.6;
        eyecam_DLC_params.lidMinStd = 7;
        eyecam_DLC_params.minDistLidCenter = 0.5; % in number of pupil heights
        eyecam_DLC_params.surroundingBlinks = 5; % in frames
        eyecam_DLC_params.interpMethod = 'linear';
        eyecam_DLC_params.maxDistPupilLid = 10; % in pixels
        eyecam_DLC_params.smoothSpan = 5; % in frames
        
        %% Some processing steps from pupil DLC page?
        
        % relation between horizontal pupil position, and width and height of pupil
        pupil_width = pupil_right(:,1) - pupil_left(:,1);
        pupil_height = pupil_bottom(:,2) - pupil_top(:,2);
        pupil_centerX = pupil_left(:,1) + 0.5 .* pupil_width;
        
        pupil_valid = all(data(:, 3:3:12) > eyecam_DLC_params.minCertainty, 2);
        
        if sum(pupil_valid) == 0
            continue
        end
        
        [F, F0] = dlc.estimateHeightFromWidthPos(pupil_width(pupil_valid), pupil_height(pupil_valid), ...
            pupil_centerX(pupil_valid), eyecam_DLC_params);
        
        [pupil_centerY_adj, pupil_height_adj] = dlc.adjustCenterHeight(data, F, eyecam_DLC_params);
        
        [blinks, bl_starts, bl_stops] = dlc.detectBlinks(data, eyecam_DLC_params, ...
            pupil_centerY_adj, pupil_height_adj, true);
        
        pupil_centerX = medfilt1(pupil_centerX, eyecam_DLC_params.smoothSpan);
        pupil_centerY_adj = medfilt1(pupil_centerY_adj, eyecam_DLC_params.smoothSpan);
        pupil_center = [pupil_centerX, pupil_centerY_adj];
        pupil_center(blinks,:) = NaN;
        pupil_height_adj = medfilt1(pupil_height_adj, eyecam_DLC_params.smoothSpan);
        diameter = pupil_height_adj;
        diameter(blinks) = NaN;
        
        % Simpler method
        
        pupil_width = eyecam_dlc.pupil_right.x - eyecam_dlc.pupil_left.x;
        pupil_height = eyecam_dlc.pupil_bot.y - eyecam_dlc.pupil_top.y;
        pupil_centerX = eyecam_dlc.pupil_left.x + 0.5 .* pupil_width;
        
        eyecam_dlc_likelihood_cell = cellfun(@(name) eyecam_dlc.(name).likelihood, fieldnames(eyecam_dlc), 'UniformOutput', false);
        eyecam_dlc_likelihood = [eyecam_dlc_likelihood_cell{:}];
        
        pupil_valid = all(eyecam_dlc_likelihood(:, :) > eyecam_DLC_params.minCertainty, 2);
        
        pupil_likelihood = eyecam_dlc_likelihood(:,5:8);
        
        uncertain_pupil_points = sum(pupil_likelihood > 0.8, 2) <= 2;        % also previous and next 5 frames
        tmp = uncertain_pupil_points;
        for t = 1:eyecam_DLC_params.surroundingBlinks
            uncertain_pupil_points = uncertain_pupil_points | [false(t,1); tmp(1:end-t)];
            uncertain_pupil_points = uncertain_pupil_points | [tmp(1+t:end); false(t,1)];
        end
        
        time_uncertain_pupil_points = eyecam_t(uncertain_pupil_points);
        time_blinks = eyecam_t(blinks);
        
        % save
        all_blinks = blinks | uncertain_pupil_points;
        both = sum(blinks & uncertain_pupil_points)/sum(all_blinks);
        mine = sum(~blinks & uncertain_pupil_points)/sum(all_blinks);
        not_mine = sum(blinks & ~uncertain_pupil_points)/sum(all_blinks);
        
        blink_info = [both mine not_mine];
        
        eyeblink(animal_id).blinks(day_index,:) = blink_info;
    end
    disp(['Done with ' animal])
end

disp('Done all')

% Save 
save('eyeblink.mat', 'eyeblink')
disp('Saved')
