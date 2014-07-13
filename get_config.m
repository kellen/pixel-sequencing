 function[config] = get_config()
    config = containers.Map;
    config('debug') = false;
    config('debug_slices') = false;
    config('debug_do_registration_final') = config('debug') || true;
    
    config('config_output_filename') = 'config.txt';

    % slide name prefix
    config('prefix') = 'slideA';
    
    % these are rough; empirically measured.
    config('position_overlap_h') = 150; % 175
    config('position_overlap_v') = 150; % 140
    config('position_row_length') = 4;
    
    config('positions') = 1:16; %5:5; %[5 12 13 14]; %1:16; %14:14;% 1:16; %4:4; %1:16;
    % problem with DO alignment on:
    %5 12 13 14
    
    config('cycles') = 1:4; % 1:4
    config('input_cycles') = {'1st', '2nd', '3rd', '4th'};
    % T=1, G=2, C=3, A=4 to correspond with ID_list_BCpanel
    % this is the order in which the image channels come, 
    % i.e. T is the first channel (channel 2), G is the 2nd (channel 3)
    config('types') = {'T', 'G', 'C', 'A'}; 
    config('channels') = 2:5;
    config('do_dirname') = 'DO';
    config('do_channel') = 2;
    
    % output filetype
    config('filetype') = 'Tiff';
    
    config('tophat_size') = 4; %10;
    
    config('do_registration_ignore_T') = true; % T channel is the worst! :(
    
    config('do_registration_tophat_size') = 4;
    
    config('do_registration_use_gaussian') = false;
    config('do_registration_gaussian_size') = 6;
    config('do_registration_gaussian_sigma') = 1;
    
    config('do_registration_use_hitandmiss') = false;
    
    config('do_registration_use_size_threshold') = true;
    config('do_registration_size_threshold') = 4;
    
    config('do_registration_cpd_opt') = struct('method','rigid','scale',0,'viz',config('debug'),'max_it',300,'outliers',0.001,'tol',1e-6);

    config('do_registration_hitandmiss_before_pointdrift_interval') = [ ...
     -1    -1    -1    -1    -1    -1    -1    -1    -1    -1    -1    -1    -1 ;
     -1     0     0     0     0     0     0     0     0     0     0     0    -1 ;
     -1     0     0     0     0     0     0     0     0     0     0     0    -1 ;
     -1     0     0     0     0     0     0     0     0     0     0     0    -1 ;
     -1     0     0     0     0     0     0     0     0     0     0     0    -1 ;
     -1     0     0     0     0     0     0     0     0     0     0     0    -1 ;
     -1     0     0     0     0     0     1     0     0     0     0     0    -1 ;
     -1     0     0     0     0     0     0     0     0     0     0     0    -1 ;
     -1     0     0     0     0     0     0     0     0     0     0     0    -1 ;
     -1     0     0     0     0     0     0     0     0     0     0     0    -1 ;
     -1     0     0     0     0     0     0     0     0     0     0     0    -1 ;
     -1     0     0     0     0     0     0     0     0     0     0     0    -1 ;
     -1    -1    -1    -1    -1    -1    -1    -1    -1    -1    -1    -1    -1 ;
     ];
    
    config('pointdrift_tophat_radius') = 1;
    % only use DoG to clean up images for pointdrift?
    config('pointdrift_use_dog') = true;
    % find pointdrift points on each base image rather than the MIP
    config('pointdrift_points_from_base_images') = false;
    % blur the images before applying the quantization?
    config('pointdrift_gaussian') = false;
    % apply quantization (false = use im2bw)
    config('pointdrift_quantize') = true;
    config('pointdrift_quantize_classes') = 3;
    
    
    
    
    % FIXME this seems sketchy
    config('registration_ignore_rotation') = false;
    
    
    % which registration type to use
    config('registration_type') = 'ransac'; % 'gradient'; %'ransac';

    
    
    
    
    config('cpd_opt') = struct('method','rigid','scale',0,'viz',config('debug'),'max_it',300,'outliers',0.001,'tol',1e-6);
    %cpd_opt = struct('method','rigid','scale',0,'viz',1,'max_it',300,'outliers',0.1,'tol',1e-5);
    
    % ransac registration configuration
    %config('surf_features_metric_threshold') = 600; % 1000.0 is default
    config('surf_features_metric_threshold') = 1000;
    config('ransac_estimation_type') = 'similarity'; % 'affine'
    
    % deprecated
    % config('register_by') = 'pair'; %'row'; % 'all'

    
    
    config('slice_registration_iterations') = 100;
    % unclear if these are relative to the pixel or the input space
    config('slice_registration_steplength_min') = .000005; % 0.00001 is default
    config('slice_registration_steplength_max') = 0.02; % 0.0625 is default
    config('slice_registration_relaxation_factor') = 0.5; % 0.5 is default
    
    
    % gradient descent registration configuration 
    % number of iterations during registration
    config('registration_iterations') = 500;
    % unclear if these are relative to the pixel or the input space
    config('registration_steplength_min') = .000005; % 0.00001 is default
    config('registration_steplength_max') = 0.06; % 0.0625 is default
    config('registration_relaxation_factor') = 0.5; % 0.5 is default
    
    % the known/valid sequences and their names
    taglist = ID_list_BCpanel();
    config('taglist') = taglist;
    % the list of acceptable sequences
    config('valid') = cell2mat(taglist(:, 2));
    
    % the chosen threshold set for this run.
    config('thresholds_key') = 'thresholds_tophat';
    
    thresholds = containers.Map;
    thresholds('avg') = 80;
    thresholds('int') = 200;
    thresholds('quality') = 0.31;
    thresholds('size') = 50;
    % thresholds for 'unmanipulated' images
    config('thresholds_raw') = thresholds;
    
    thresholds_tophat = containers.Map;
    thresholds_tophat('avg') = 25; %10;
    thresholds_tophat('int') = 40; %100;
    thresholds_tophat('quality') = 0.475; % 0.52;
    thresholds_tophat('size_lower') = 5;
    thresholds_tophat('size_upper') = 80; % 50;
    % thresholds for tophat-transformed images
    config('thresholds_tophat') = thresholds_tophat;
    
    % input dir/file formats
    config('dir_input') = '%s_%d_%s';
    config('img_infile') = [config('dir_input') '%s_%d_%s_c%d.tif'];
    
    % output
    config('dir_prefix') = config('prefix');
    config('dir_position') = [config('dir_prefix') '/%d'];
    config('dir_cycle') = [config('dir_position') '/%d'];
    config('img_cycle') = [config('dir_cycle') '/%s.tif'];
    config('img_do') = [config('dir_position') '/do.tif'];
    config('img_mip') = [config('dir_position') '/%d-mip.tif'];
    config('img_edge') = [config('dir_position') '/edge-%d.tif'];
    config('img_dist') = [config('dir_position') '/dist-%d.tif'];
    config('img_do_dist') = [config('dir_position') '/dist-do.tif'];
    
    config('mat_sequence') = [config('dir_position') '/sequence.mat'];
    % the current threshold map to appropriately name the output file
    t = config(config('thresholds_key'));
    config('img_filtered') = sprintf('%s/%s', config('dir_position'), ...
        sprintf('filtered_avg-%d-int-%d-quality-%f-lsize-%d-usize-%d.tif', ...
        t('avg'), t('int'), t('quality'), t('size_lower'), t('size_upper')));
    
    config('mat_centroids') = [config('dir_position') '/centroids.mat'];
    config('mat_cellprofiler') = 'DefaultOUT.mat';
    
    % label format for cellprofiler output matrix
    config('label_cellprofiler') = 'Intensity_MaxIntensity_%d_%d';
    
    % quality threshold for the cellprofiler blobs
    config('threshold_quality_cellprofiler') = 0.5;
    
    
    % quality steps
    config('show_steps_detail_range') = 0:300;
    config('show_steps_quality') = 0:0.1:1.0;
    config('show_steps_size') = 15:30; % 0:15
    config('show_steps_size_upper') = 100:10:150; %10:10:100;
    config('show_steps_average') = 80:4:120; %0:4:80;
    config('show_steps_intensity') = 0:10:160;
    
    % ROI
    config('roi_x') = 500:800;
    config('roi_y') = 250:500;
    
    % don't compare these
    config('skip') = [5 12 13 14];
end