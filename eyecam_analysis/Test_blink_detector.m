eyecam_dlc_params.surroundingBlinks = 5;
eyecam_dlc_params.thresLikelihood = 0.8;

pupil_likelihood = [ ...
    eyecam_dlc.pupil_bot.likelihood, ...
    eyecam_dlc.pupil_left.likelihood, ...
    eyecam_dlc.pupil_right.likelihood, ...
    eyecam_dlc.pupil_top.likelihood];

uncertain_pupil_points = conv(sum(pupil_likelihood > ...
    eyecam_dlc_params.thresLikelihood, 2) <= 2, ... 
    ones(2*eyecam_dlc_params.surroundingBlinks+1,1), 'same') > 0;

pupil_diameter = ...
    (vecnorm([eyecam_dlc.pupil_right.x - eyecam_dlc.pupil_left.x, ...
    eyecam_dlc.pupil_right.y - eyecam_dlc.pupil_left.y]') + ...
    vecnorm([eyecam_dlc.pupil_top.x - eyecam_dlc.pupil_bot.x, ...
    eyecam_dlc.pupil_top.y - eyecam_dlc.pupil_bot.y]'))/2; 

pupil_diameter(uncertain_pupil_points) = NaN;
