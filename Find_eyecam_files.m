
%% Paths to stuff
addpath(genpath(cd));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\npy-matlab'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\widefield'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\Lilrig'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\PupilDetection_DLC'));
addpath(genpath('C:\Users\Andrada\Documents\GitHub\AP_scripts_cortexlab'));

%% code
animal = 'AP108';
protocol = 'AP_lcrGratingPassive';
experiments = AP_find_experiments(animal, protocol);

eyecam_all_fn = {};
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

formatted_eyecam_all_fn = cellfun(@(x) ['r''' x ''''], eyecam_all_fn, 'Uni', 0);

fid = fopen('eyecam_paths.txt','w');

fprintf(fid, 'videos = [%s]', strjoin(formatted_eyecam_all_fn, ', '))

fclose(fid)