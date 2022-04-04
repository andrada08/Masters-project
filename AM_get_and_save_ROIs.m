%% Get and save ROIs

roi = struct;

animal = 'AP107';
day = '2021-12-16';
experiment = 1;
verbose = true;
AP_load_experiment;

% temp
Udf_aligned = AP_align_widefield(Udf,animal,day);

% frontal left
avg_im_aligned = AP_align_widefield(avg_im,animal,day);
figure;
imagesc(avg_im_aligned);
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
axis image off

this_roi = roipoly();
roi.frontal_left = this_roi;


% frontal right
avg_im_aligned = AP_align_widefield(avg_im,animal,day);
figure;
imagesc(avg_im_aligned);
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
axis image off

this_roi = roipoly();
roi.frontal_right = this_roi;

% temp - so confused?????
left_trace = AP_svd_roi(Udf_aligned,fVdf,[],[],roi.frontal_left);
right_trace = AP_svd_roi(Udf_aligned,fVdf,[],[],roi.frontal_right);

figure;
plot(left_trace, 'b')
hold on
plot(right_trace, 'r-')

% visual left
avg_im_aligned = AP_align_widefield(avg_im,animal,day);
figure;
imagesc(avg_im_aligned);
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
axis image off

this_roi = roipoly();

roi.visual_left = this_roi;

% visual right
avg_im_aligned = AP_align_widefield(avg_im,animal,day);
figure;
imagesc(avg_im_aligned);
AP_reference_outline('ccf_aligned',[0.5,0.5,0.5]);
axis image off

this_roi = roipoly();

roi.visual_right = this_roi;

save('ROIs.mat', 'roi');