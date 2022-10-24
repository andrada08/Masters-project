
%% Paths to stuff
addpath(genpath(cd));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\npy-matlab'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\widefield'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\Lilrig'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\PupilDetection_DLC'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\AP_scripts_cortexlab'));

%% code
% animals = {'AP100','AP101','AP103','AP104','AP105','AP106','AP107','AP108','AP109'};

animals = {'AP107','AP108', 'AP109', 'AP110', 'AP111', 'AP112', 'AP113', 'AP114', 'AP115'};

eyecam_all_fn = {};
for animal_id=1:length(animals)
    animal = animals{animal_id};
    experiments = AP_find_experiments(animal);
    for day_index=1:length(experiments)
        day = experiments(day_index).day;
        exp_day = AP_list_experiments(animal,day);
        for experiment=[exp_day.experiment]
            [eyecam_fn, eyecam_exists] = AP_cortexlab_filename(animal,day,experiment,'eyecam');
            if ~eyecam_exists
                continue
            end
            eyecam_all_fn = [eyecam_all_fn {eyecam_fn}];
        end
    end
end

formatted_eyecam_all_fn = cellfun(@(x) ['r''' x ''''], eyecam_all_fn, 'Uni', 0);

fid = fopen('eyecam_paths.txt','w');

fprintf(fid, 'videos = [%s]', strjoin(formatted_eyecam_all_fn, ', '))

fclose(fid)

%% Find and load early and late video for each mouse

animals = {'AP100','AP101','AP103','AP104','AP105','AP106','AP107','AP108','AP109'};

for animal_id=1:length(animals)
    animal = animals{animal_id};
    experiments = AP_find_experiments(animal);
    
    % early day - day 4
    day = experiments(4).day;
    experiment = 1;
    load_parts.imaging = false;
    load_parts.cam = true;
    AP_load_experiment;
    
    AP_mousemovie(eyecam_fn,eyecam_dlc)
    
    % late day - second to last day
    day = experiments(length(experiments)-2).day;
    experiment = 1;
    
    % there is just one session from AP101 that doesn't have eyecam 
    [eyecam_fn, eyecam_exists] = AP_cortexlab_filename(animal,day,experiment,'eyecam');
    if ~eyecam_exists
        day = experiments(length(experiments)-3).day;
    end
    
    load_parts.imaging = false;
    load_parts.cam = true;
    AP_load_experiment;
    
    AP_mousemovie(eyecam_fn,eyecam_dlc) 
    
end