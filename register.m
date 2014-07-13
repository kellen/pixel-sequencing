function [] = register(input_dir, mip_input_dir, output_dir) 
    config = get_config();
    write_config(output_dir);

    regtype = config('registration_type');
    disp(['Registering with ' regtype '...']);
    if strcmp(regtype, 'ransac')
        ransac(input_dir, mip_input_dir, output_dir);
    elseif strcmp(regtype, 'gradient') || strcmp(regtype, 'corr')
        gradient_descent(input_dir, mip_input_dir, output_dir);
    elseif strcmp(regtype, 'bwmutual')
        bwmutual(input_dir, mip_input_dir, output_dir);
    elseif strcmp(regtype, 'bwransac')
        bwransac(input_dir, mip_input_dir, output_dir);
    elseif strcmp(regtype, 'bwcorr')
        bwcorr(input_dir, mip_input_dir, output_dir);
    elseif strcmp(regtype, 'distreg')
        distreg(input_dir, mip_input_dir, output_dir);
    elseif strcmp(regtype, 'ptsransac')
        ptsransac(input_dir, mip_input_dir, output_dir);
    elseif strcmp(regtype, 'pointdrift')
        pointdrift(input_dir, mip_input_dir, output_dir);
    end
    
end

function [] = pointdrift(input_dir, mip_input_dir, output_dir)
    config = get_config();
    positions = config('positions');
    cpd_opt = config('cpd_opt');
    
    for posidx = 1:numel(positions) 
        disp(['Position ' num2str(posidx) '...']);
        position = positions(posidx);
        
        pos = read_position(input_dir, position);
        c = pos('cycles');
        m = read_mips(mip_input_dir, position);
        
        % register against the DO
        do_reg = pos('do');
        cents = pointdrift_points(do_reg);
        sz = imref2d(size(do_reg));
        
        registered = cell([numel(c) 1]);
        registered_mips = cell([numel(c) 1]);
        for cycleidx = 1:numel(c)
            %mcents = pointdrift_points(m{cycleidx});
            %[trans,~] = cpd_register(cents, mcents, cpd_opt);
            %tform = affine2d([trans.R [0;0]; [trans.t' 1]]);
            if config('pointdrift_points_from_base_images')
                zsize = size(c{cycleidx},3);
                mcents = [];
                cur = c{cycleidx};
                for baseidx=1:zsize
                    mcents = [mcents ; pointdrift_points(cur(:,:,baseidx))];
                end
            else 
                mcents = pointdrift_points(m{cycleidx});
            end
            
            if config('debug') ; 
                figure;
                imshow(m{cycleidx});
                hold on,scatter(cents(:, 1), cents(:, 2), 15, repmat([0,0,1] ,size(cents,1),1), '+', 'red'); 
                hold on,scatter(mcents(:, 1), mcents(:, 2), 15, repmat([0,0,1] ,size(mcents,1),1), 'o', 'blue'); 
                title(['Centroids ' num2str(cycleidx)]);
            end
            [trans,~] = cpd_register(cents, mcents, cpd_opt);
            tform = affine2d([trans.R [0;0]; [trans.t' 1]]);
            
            %keyboard;
            %figure,imshow(do_reg),hold on,cpd_plot_iter(cents, mcents); title(['Before: ' num2str(cycleidx)]);
            %figure,imshow(do_reg),hold on,cpd_plot_iter(cents, trans.Y);  title(['After: '  num2str(cycleidx)]);
            
            % register all the channels for this cycle at once
            % imwarp results in size(fixed) images, but all the inputs should be the exact same size
            registered{cycleidx} = imwarp(c{cycleidx}, tform, 'cubic', 'OutputView', sz);
            % register the mips
            registered_mips{cycleidx} = imwarp(m{cycleidx}, tform, 'cubic', 'OutputView', sz);
            if config('debug')
                figure; imshow(imfuse(do_reg, registered_mips{cycleidx}));
                title(['Registered ' num2str(cycleidx)]);
            end
            
            % clamp values to strictly positive
            registered{cycleidx}(registered{cycleidx} < 0) = 0;
            
            disp(['Registered cycle ' num2str(cycleidx) '.']);
        end
        save_images(position, c, registered, registered_mips, do_reg, output_dir)
    end
end

function [] = ptsransac(input_dir, mip_input_dir, output_dir)
    config = get_config();
    positions = config('positions');
    for posidx = 1:numel(positions) 
        disp(['Position ' num2str(posidx) '...']);
        position = positions(posidx);
        
        pos = read_position(input_dir, position);
        c = pos('cycles');
        m = read_mips(mip_input_dir, position);
        
        % register against the DO
        do_reg = pos('do');
        % use features from the first mip
        bw = im2bw(pos('do'));
        points = detectSURFFeatures(bw, 'MetricThreshold', config('surf_features_metric_threshold'));
        [features, valid] = extractFeatures(bw,  points);
        sz = imref2d(size(bw));
        
        registered = cell([numel(c) 1]);
        registered_mips = cell([numel(c) 1]);
        for cycleidx = 1:numel(c)            
            r_points = detectSURFFeatures(im2bw(m{cycleidx}), 'MetricThreshold', config('surf_features_metric_threshold'));
            [r_features, r_valid] = extractFeatures(im2bw(m{cycleidx}), r_points);
            pairs = matchFeatures(features, r_features);

            matched  = valid(pairs(:,1));
            r_matched = r_valid(pairs(:,2));
            [tform,~,~] = estimateGeometricTransform(r_matched,matched,config('ransac_estimation_type'));

            % register all the channels for this cycle at once
            % imwarp results in size(fixed) images, but all the inputs should be the exact same size
            registered{cycleidx} = imwarp(c{cycleidx}, tform, 'cubic', 'OutputView', sz);
            % register the mips
            registered_mips{cycleidx} = imwarp(m{cycleidx}, tform, 'cubic', 'OutputView', sz);
            % clamp values to strictly positive
            registered{cycleidx}(registered{cycleidx} < 0) = 0;
            
            disp(['Registered cycle ' num2str(cycleidx) '.']);
        end
        save_images(position, c, registered, registered_mips, do_reg, output_dir)
    end
end

function [] = distreg(input_dir, mip_input_dir, output_dir)
    config = get_config();
    cycles = config('cycles');
    
        % images produced on same device
        [optimizer, metric] = imregconfig('monomodal');
        optimizer.MaximumIterations = config('registration_iterations');
        optimizer.MaximumStepLength = config('registration_steplength_max');
        optimizer.MinimumStepLength = config('registration_steplength_min');
        optimizer.RelaxationFactor = config('registration_relaxation_factor');
        disp(['Registering images with ' num2str(optimizer.MaximumIterations) ' iterations']);

    positions = config('positions');
    for posidx = 1:numel(positions) 
        disp(['Position ' num2str(posidx) '...']);
        position = positions(posidx);
        
        pos = read_position(input_dir, position);
        c = pos('cycles');
        m = read_mips(mip_input_dir, position);
        
        do_reg = pos('do');
        fixed = imcomplement(bwdist(im2bw(do_reg)));
        
        % use the first mip as the base image
        registered = cell([numel(c) 1]);
        registered_mips = cell([numel(c) 1]);
        dists = cell([numel(c) 1]);
        for cycleidx = 1:numel(c)
            moving = imcomplement(bwdist(im2bw(m{cycleidx})));

            % should only be rigid transformations i.e. translation and rotation
            tform = imregtform(moving, fixed, 'rigid', optimizer, metric);
            dists{cycleidx} = moving;

            % register all the channels for this cycle at once
            % imwarp results in size(fixed) images, but all the inputs should be the exact same size
            registered{cycleidx} = imwarp(c{cycleidx}, tform, 'cubic', 'OutputView', imref2d(size(fixed)));
            % register the mips
            registered_mips{cycleidx} = imwarp(m{cycleidx}, tform, 'cubic', 'OutputView', imref2d(size(fixed)));
            % clamp values to strictly positive
            registered{cycleidx}(registered{cycleidx} < 0) = 0;
            disp(['Registered cycle ' num2str(cycleidx) '.']);
        end
        save_images(position, c, registered, registered_mips, do_reg, output_dir);
        fmt = '%s/%s';
        fn = sprintf(fmt, output_dir, sprintf(config('img_do_dist'), position));
        write_image(double(fixed), fn);
        for cycleidx = 1:numel(c)
            cycle = cycles(cycleidx);
            % write the registered edge images
            lr = dists{cycleidx};
            fn = sprintf(fmt, output_dir, sprintf(config('img_dist'), position, cycle));
            write_image(double(lr), fn);
        end
    end
end

function [] = bwransac(input_dir, mip_input_dir, output_dir)
    config = get_config();
    positions = config('positions');
    for posidx = 1:numel(positions) 
        disp(['Position ' num2str(posidx) '...']);
        position = positions(posidx);
        
        pos = read_position(input_dir, position);
        c = pos('cycles');
        m = read_mips(mip_input_dir, position);
        
        % register against the DO
        do_reg = pos('do');
        % use features from the first mip
        bw = im2bw(pos('do'));
        points = detectSURFFeatures(bw, 'MetricThreshold', config('surf_features_metric_threshold'));
        [features, valid] = extractFeatures(bw,  points);
        sz = imref2d(size(bw));
        
        registered = cell([numel(c) 1]);
        registered_mips = cell([numel(c) 1]);
        for cycleidx = 1:numel(c)            
            r_points = detectSURFFeatures(im2bw(m{cycleidx}), 'MetricThreshold', config('surf_features_metric_threshold'));
            [r_features, r_valid] = extractFeatures(im2bw(m{cycleidx}), r_points);
            pairs = matchFeatures(features, r_features);

            matched  = valid(pairs(:,1));
            r_matched = r_valid(pairs(:,2));
            [tform,~,~] = estimateGeometricTransform(r_matched,matched,'affine');

            % register all the channels for this cycle at once
            % imwarp results in size(fixed) images, but all the inputs should be the exact same size
            registered{cycleidx} = imwarp(c{cycleidx}, tform, 'cubic', 'OutputView', sz);
            % register the mips
            registered_mips{cycleidx} = imwarp(m{cycleidx}, tform, 'cubic', 'OutputView', sz);
            % clamp values to strictly positive
            registered{cycleidx}(registered{cycleidx} < 0) = 0;
            
            disp(['Registered cycle ' num2str(cycleidx) '.']);
        end
        save_images(position, c, registered, registered_mips, do_reg, output_dir)
    end
end

function [] = ransac(input_dir, mip_input_dir, output_dir)
    config = get_config();
    positions = config('positions');
    for posidx = 1:numel(positions) 
        disp(['Position ' num2str(posidx) '...']);
        position = positions(posidx);
        
        pos = read_position(input_dir, position);
        c = pos('cycles');
        m = read_mips(mip_input_dir, position);
        
        % use the first mip as the base image
        registered = cell([numel(c) 1]);
        registered{1} = c{1};
        registered_mips = cell([numel(c) 1]);
        registered_mips{1} = m{1};
        
        % use features from the first mip
        points = detectSURFFeatures(m{1}, 'MetricThreshold', config('surf_features_metric_threshold'));
        [features, valid] = extractFeatures(m{1},  points);
        sz = imref2d(size(m{1}));
        
        for cycleidx = 2:numel(c)
            tform = ransac_reg(features, valid, m{cycleidx});

            % register all the channels for this cycle at once
            % imwarp results in size(fixed) images, but all the inputs should be the exact same size
            registered{cycleidx} = imwarp(c{cycleidx}, tform, 'cubic', 'OutputView', sz);
            % register the mips
            registered_mips{cycleidx} = imwarp(m{cycleidx}, tform, 'cubic', 'OutputView', sz);
            % clamp values to strictly positive
            registered{cycleidx}(registered{cycleidx} < 0) = 0;
            disp(['Registered cycle ' num2str(cycleidx) '.']);
        end

        [do_reg] = register_do(pos('do'), registered);
            
        save_images(position, c, registered, registered_mips, do_reg, output_dir)
    end
end

function [do_reg] = register_do(do, registered)
    config = get_config();
    mip_mip = do_registration_mip(registered);
    if config('debug')
        figure,imshow(mip_mip);title('MIP of MIPs'); 
    end
    
    if config('do_registration_use_gaussian')
        h = fspecial('gaussian', config('do_registration_gaussian_size'), config('do_registration_gaussian_sigma'));
        mip_mip = imfilter(mip_mip, h);
        if config('debug')
            figure,imshow(mip_mip);title('Gaussian of MIP of MIPs'); 
        end
    end
    
    mip_mip_th = imtophat(mip_mip, strel('disk', config('do_registration_tophat_size')));
    do_th = imtophat(do, strel('disk', config('do_registration_tophat_size')));
    
    if config('debug')
        %figure,imshow(mip_mip_th);title('Tophat of MIP of MIPs'); 
        %figure,imshow(do_th);title('Tophat of DO'); 
        figure,imshow(imfuse(do_th, mip_mip_th));title('TH of DO with TH of MIP of MIPS');
    end
    
    tform = do_pointdrift_reg(do_th, mip_mip_th);
    do_reg = imwarp(do, tform, 'cubic', 'OutputView', imref2d(size(do)));
    do_reg(do_reg < 0) = 0;
    
    if config('debug_do_registration_final')
        figure,imshow(imfuse(mip_mip, do_reg));title('Aligned: DO + Tophat of MIP of MIPs'); 
    end
end

function [mippy] = do_registration_mip(registered)
    config = get_config();
    cycles = config('cycles');
    types = config('types');
    
    mippy = uint8(zeros(size(registered{1}(:,:,1))));
    
    for cycleidx = 1:numel(cycles)
        if config('do_registration_ignore_T')
            for typeidx = 1:numel(types)
                if ~strcmp(types{typeidx}, 'T')
                    mippy = max(mippy, registered{cycleidx}(:,:,typeidx));
                end
            end
        else 
            mippy = max(mippy, max(registered{cycleidx}, [], 3));
        end
    end
end

function [tform] = do_pointdrift_reg(im, im_reg)
    config = get_config();
    cpd_opt = config('do_registration_cpd_opt');
    
    impts = do_pointdrift_reg_points(im);
    regpts = do_pointdrift_reg_points(im_reg);
    
    if config('debug') ; 
        figure;
        imshow(im_reg);
        hold on,scatter(regpts(:, 1), regpts(:, 2), 15, repmat([0,0,1] ,size(regpts,1),1), '+', 'red'); 
        hold on,scatter(impts(:, 1), impts(:, 2), 15, repmat([0,0,1] ,size(impts,1),1), 'o', 'blue'); 
        title('Centroids pre-registration');
    end
    [trans,~] = cpd_register(regpts, impts, cpd_opt);
    tform = affine2d([trans.R [0;0]; [trans.t' 1]]);
    if config('debug') ; 
        figure;
        imshow(im_reg);
        hold on,scatter(regpts(:, 1), regpts(:, 2), 15, repmat([0,0,1] ,size(regpts,1),1), '+', 'red'); 
        hold on,scatter(trans.Y(:, 1), trans.Y(:, 2), 15, repmat([0,0,1] ,size(trans.Y,1),1), 'o', 'blue'); 
        title('Centroids post-registration');
    end
end

function [cents] = do_pointdrift_reg_points(im)
    config = get_config();
    
    % quantize into 3 classes
    numclasses = config('pointdrift_quantize_classes');
    quant = imquantize(im, multithresh(im, numclasses-1));
    % grab class 3, which is (presumably) the brightest points
    bw = (quant==numclasses);
    t = ['quant==' num2str(numclasses)];
    
    if config('do_registration_use_hitandmiss')
       bw = bw & bwhitmiss(quant>1, config('do_registration_hitandmiss_before_pointdrift_interval'));
       t = [t ' & hist&miss'];
    end
    
    if config('do_registration_use_size_threshold')
        bw = bwareaopen(bw, config('do_registration_size_threshold'));
        %bw(rm) = 0;
        t = [t ' & rm < ' num2str(config('do_registration_size_threshold'))];
    end
    
    if config('debug') ; 
        figure,imshow(label2rgb(quant, 'lines', 'k'));title('quant'); 
        figure,imshow(bw);title(t);
    end
    
    props = regionprops(bw, 'Centroid');
    cents = cat(1, props.Centroid);
    cents = cents(isfinite(cents(:, 1)), :);
    if config('debug')
        figure,imshow(im),hold on,scatter(cents(:, 1), cents(:, 2), 15, repmat([0,0,1] ,size(cents,1),1), '+', 'blue');title('points');
    end
end

function [tform] = ransac_reg(features, valid, img) 
    config = get_config();
    r_points = detectSURFFeatures(img, 'MetricThreshold', config('surf_features_metric_threshold'));
    [r_features, r_valid] = extractFeatures(img, r_points);
    pairs = matchFeatures(features, r_features);

    matched  = valid(pairs(:,1));
    r_matched = r_valid(pairs(:,2));
    [tform,~,~] = estimateGeometricTransform(r_matched,matched,config('ransac_estimation_type'));
end

function [] = bwmutual(input_dir, mip_input_dir, output_dir)
    config = get_config();
    cycles = config('cycles');
    
    [optimizer, metric] = imregconfig('multimodal');
    
    positions = config('positions');
    for posidx = 1:numel(positions) 
        disp(['Position ' num2str(posidx) '...']);
        position = positions(posidx);
        
        pos = read_position(input_dir, position);
        c = pos('cycles');
        m = read_mips(mip_input_dir, position);
        
        % use the first mip as the base image
        registered = cell([numel(c) 1]);
        registered{1} = c{1};
        registered_mips = cell([numel(c) 1]);
        registered_mips{1} = m{1};
        
            % register against the laplacian
            fixed = im2bw(m{1});
            registered_laps = cell([numel(c) 1]);
            registered_laps{1} = fixed;
        
        for cycleidx = 2:numel(c)
                moving = im2bw(m{cycleidx});
                % should only be rigid transformations i.e. translation and rotation
                tform = imregtform(moving, fixed, 'rigid', optimizer, metric);
            
            % register all the channels for this cycle at once
            % imwarp results in size(fixed) images, but all the inputs should be the exact same size
            registered{cycleidx} = imwarp(c{cycleidx}, tform, 'cubic', 'OutputView', imref2d(size(fixed)));
            % register the mips
            registered_mips{cycleidx} = imwarp(m{cycleidx}, tform, 'cubic', 'OutputView', imref2d(size(fixed)));
                % register the laps
                registered_laps{cycleidx} = imwarp(im2bw(m{cycleidx}), tform, 'cubic', 'OutputView', imref2d(size(fixed)));
            % clamp values to strictly positive
            registered{cycleidx}(registered{cycleidx} < 0) = 0;
            
            disp(['Registered cycle ' num2str(cycleidx) '.']);
        end
        
        % register the DO
            moving = im2bw(pos('do'));
            tform = imregtform(moving, fixed, 'rigid', optimizer, metric);

        do_reg = imwarp(pos('do'), tform, 'cubic', 'OutputView', imref2d(size(fixed)));
        do_reg(do_reg < 0) = 0;
        
        save_images(position, c, registered, registered_mips, do_reg, output_dir);
        
            fmt = '%s/%s';
            for cycleidx = 1:numel(c)
                cycle = cycles(cycleidx);
                % write the registered edge images
                lr = registered_laps{cycleidx};
                fn = sprintf(fmt, output_dir, sprintf(config('img_edge'), position, cycle));
                write_image(lr, fn);
            end
    end
end

function [] = bwcorr(input_dir, mip_input_dir, output_dir)
    config = get_config();
    cycles = config('cycles');

    positions = config('positions');
    for posidx = 1:numel(positions) 
        disp(['Position ' num2str(posidx) '...']);
        position = positions(posidx);
        
        pos = read_position(input_dir, position);
        c = pos('cycles');
        m = read_mips(mip_input_dir, position);
        
        % register the DO
        do_reg = pos('do');
        fixed = im2bw(do_reg);
        
        % use the first mip as the base image
        registered = cell([numel(c) 1]);
        registered_mips = cell([numel(c) 1]);

        for cycleidx = 1:numel(c)
            moving = im2bw(m{cycleidx});
            tform = imregcorr(moving, fixed, 'rigid');
            
            % register all the channels for this cycle at once
            % imwarp results in size(fixed) images, but all the inputs should be the exact same size
            registered{cycleidx} = imwarp(c{cycleidx}, tform, 'cubic', 'OutputView', imref2d(size(fixed)));
            % register the mips
            registered_mips{cycleidx} = imwarp(m{cycleidx}, tform, 'cubic', 'OutputView', imref2d(size(fixed)));
            % clamp values to strictly positive
            registered{cycleidx}(registered{cycleidx} < 0) = 0;
            
            disp(['Registered cycle ' num2str(cycleidx) '.']);
        end
        save_images(position, c, registered, registered_mips, do_reg, output_dir);
    end
end

function [] = gradient_descent(input_dir, mip_input_dir, output_dir)
    config = get_config();
    
    % images produced on same device
    [optimizer, metric] = imregconfig('monomodal');
    optimizer.MaximumIterations = config('registration_iterations');
    optimizer.MaximumStepLength = config('registration_steplength_max');
    optimizer.MinimumStepLength = config('registration_steplength_min');
    optimizer.RelaxationFactor = config('registration_relaxation_factor');
    disp(['Registering images with ' num2str(optimizer.MaximumIterations) ' iterations']);
    
    positions = config('positions');
    for posidx = 1:numel(positions) 
        disp(['Position ' num2str(posidx) '...']);
        position = positions(posidx);
        
        pos = read_position(input_dir, position);
        c = pos('cycles');
        m = read_mips(mip_input_dir, position);
           
        % 
        do_reg = pos('do');
        fixed = do_reg;
        
        registered = cell([numel(c) 1]);
        registered_mips = cell([numel(c) 1]);
        for cycleidx = 1:numel(c)
            moving = m{cycleidx};
            % should only be rigid transformations i.e. translation and rotation
            tform = imregtform(moving, fixed, 'rigid', optimizer, metric);
            
            % register all the channels for this cycle at once
            % imwarp results in size(fixed) images, but all the inputs should be the exact same size
            registered{cycleidx} = imwarp(c{cycleidx}, tform, 'cubic', 'OutputView', imref2d(size(fixed)));
            % register the mips
            registered_mips{cycleidx} = imwarp(m{cycleidx}, tform, 'cubic', 'OutputView', imref2d(size(fixed)));
            % clamp values to strictly positive
            registered{cycleidx}(registered{cycleidx} < 0) = 0;
            
            disp(['Registered cycle ' num2str(cycleidx) '.']);
        end
        save_images(position, c, registered, registered_mips, do_reg, output_dir);
    end
end

function [] = save_images(position, c, registered, registered_mips, do_reg, output_dir) 
    config = get_config();
    cycles = config('cycles');
    types = config('types');

    % write out all images
    fmt = '%s/%s';
    for cycleidx = 1:numel(c)
        cycle = cycles(cycleidx);
        % write the registered images
        r = registered{cycleidx};
        for channel = 1:size(r, 3)
            fn = sprintf(fmt, output_dir, sprintf(config('img_cycle'), position, cycle, types{channel}));
            write_image(r(:,:,channel), fn);
        end
        % write the registred mips
        mr = registered_mips{cycleidx};
        fn = sprintf(fmt, output_dir, sprintf(config('img_mip'), position, cycle));
        write_image(mr, fn);

    end
    % write the do
    fn = sprintf(fmt, output_dir, sprintf(config('img_do'), position));
    write_image(do_reg, fn);
end

%{
im12_1 = imread('06-mips/slideA/12/1-mip.tif');
im12_do = imread('05-tophats/slideA/12/do.tif');
im = im12_1;
do = im12_do;
dobw = im2bw(do);
doL = bwlabel(dobw);
imbw = im2bw(im);
imL = bwlabel(imbw);
doProps = regionprops(doL, 'Centroid', 'Area');
imProps = regionprops(imL, 'Centroid', 'Area');
docents = cat(1, doProps.Centroid);
imcents = cat(1, imProps.Centroid);


% swedish version

[R,T] = icp(docents, imcents);
RR = repmat(R(1,:),size(imcents,1),1);
newim = (RR .* imcents) + repmat(T',size(imcents,1),1);

%???   newim = R*imcents + T;


sz = 20;
blue = [0,0,1];
mag = [1,0,1];
red = [1,0,0];
figure, imshow(do),
hold on, scatter(docents(:, 1), docents(:, 2), sz, repmat(blue,size(docents,1),1), '+'),
hold on, scatter(imcents(:, 1), imcents(:, 2), sz, repmat(mag,size(imcents,1),1), '+'),
hold on, scatter(newim(:, 1), newim(:, 2), sz, repmat(red,size(newim,1),1), 'o');



% danish version
docents = cat(1, doProps.Centroid);
imcents = cat(1, imProps.Centroid);
docents(:, 3) = 1;
imcents(:, 3) = 1;
[R,T] = icp(docents', imcents',100,'Matching','kDtree');

newim = R*(imcents') + repmat(T,1,length(imcents'));
newim = newim';
%}
