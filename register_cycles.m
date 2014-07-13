function [] = register_cycles(input_dir, mip_input_dir, output_dir)
    config = get_config();
    write_config(output_dir);
    positions = config('positions');
    for posidx = 1:numel(positions)
        disp(['Position ' num2str(posidx) '...']);
        position = positions(posidx);
        register_position(position, input_dir, mip_input_dir, output_dir);
    end
end

function [slice,slicewith,to_slice] = slice_position(pos, withpos, expanded, withexpanded, imx, imy)
    if pos(1) == withpos(1)
        % same y coordinate == same row
        if pos(2) > withpos(2)
            % left side of pos, right side of withpos
            [slice,~] = left(pos, expanded, imx, imy);
            [slicewith,to_slice] = right(withpos, withexpanded, imx, imy);
        else
            % right side of pos, left side of withpos
            [slice,~] = right(pos, expanded, imx, imy);
            [slicewith,to_slice] = left(withpos, withexpanded, imx, imy);
        end
    else
        if pos(2) > withpos(2)
            % top of pos, bottom of withpos
            [slice,~] = top(pos, expanded, imx, imy);
            [slicewith,to_slice] = bottom(withpos, withexpanded, imx, imy);
        else
            % bottom of pos, top of withpos
            [slice,~] = bottom(pos, expanded, imx, imy);
            [slicewith,to_slice] = top(withpos, withexpanded, imx, imy);
        end
    end
    config = get_config();
    if config('debug_slices') ; figure,imshow(imfuse(slice, slicewith, 'montage')); title('slices'); end
end

function [to] = trans(sy, sx)
    to = [1 0 0; 0 1 0; sy sx 1];
end

function [from] = itrans(to)
    % assumes no rotation
    from = [1 0 0; 0 1 0; -to(3,1) -to(3,2) 1];
end

function [slice, to_slice] = bottom(pos, expanded, imx, imy)
    config = get_config();
    ov = config('position_overlap_v');

    sy = ((pos(1) + 1) * imy)+1-((pos(1)+1) * ov);
    sx = (pos(2) * imx)+1;
    to_slice = trans(sy, sx);
    slice = expanded(...
                    sy:sy+ov-1, ...
                    sx:sx+imx-1 ...
                    );
    if config('debug_slices') ; figure, imshow(slice); title('bottom slice'); end  
end

function [slice, to_slice] = top(pos, expanded, imx, imy)
    config = get_config();
    ov = config('position_overlap_v');

    sy = (pos(1) * imy)+1-(pos(1) * ov);
    sx = (pos(2) * imx)+1;
    to_slice = trans(sy, sx);
    slice = expanded( ...
                    sy:sy+ov-1, ...
                    sx:sx+imx-1 ...
                    );
    if config('debug_slices') ; figure, imshow(slice); title('top slice'); end                
end

function [slice, to_slice] = right(pos, expanded, imx, imy)
    config = get_config();
    oh = config('position_overlap_h');

    sy = (pos(1) * imy)+1;
    sx = ((pos(2) + 1) * imx)+1-((pos(2)+1) * oh);
    to_slice = trans(sy,sx);
    slice = expanded( ...
                    sy:sy+imy-1, ...
                    sx:sx+oh-1 ...
                    );
    if config('debug_slices') ; figure, imshow(slice); title('right slice'); end
end

function [slice, to_slice] = left(pos, expanded, imx, imy) 
    config = get_config();
    oh = config('position_overlap_h');

    sy = (pos(1) * imy)+1;
    sx = (pos(2) * imx)+1-(pos(2) * oh);
    to_slice = trans(sy,sx);

    slice = expanded(...
                    sy:sy + imy - 1, ...
                    sx:sx + oh - 1 ...
                    );
    if config('debug_slices') ; figure, imshow(slice); title('left slice'); end
end

function [topos,frompos] = trans_to_position(pos, imx, imy)
    config = get_config();

    ypos = pos(1);
    xpos = pos(2);

    ov = config('position_overlap_v');
    oh = config('position_overlap_h');

    topos = [1 0 0; ...
             0 1 0; ...
             ((xpos * imx) - (xpos * oh)) ((ypos * imy) - (ypos * ov)) 1];
    frompos = [1 0 0; 0 1 0; -topos(3,1) -topos(3,2) 1];
end

function [refsize,prev,pos,post] = position_in_combined(position, imx, imy)
    % pos = [y,x]
    [~, ~, oddrow, first_in_row, ~, ~, last_in_row] = position_info(position);    

    % normal case: middle of 3 items
    refsize = imref2d([imy (3*imx)]);
    prev = [0,0];
    pos = [0,1];
    post = [0,1];

    % corner cases (ha ha)
    if first_in_row || last_in_row
        % 2x2
        refsize = imref2d([(2*imy) (2*imx)]);
        if position == 1
            refsize = imref2d([imy (2*imx)]);
            prev = NaN;
            pos = [0,0];
            post = [0,1];
        elseif position == 16
            refsize = imref2d([imy (2*imx)]);
            prev = [0,0];
            pos = [0,1];
            post = NaN;
        else
            prev = [0,0];
            pos = [1,0];
            post = [1,1];
            if last_in_row
                prev = [0,0];
                pos = [0,1];
                post = [1,1];
            end
            if oddrow
                prev = [0,1];
                pos = [1,1];
                post = [1,0];
                if last_in_row
                    prev = [0,1];
                    pos = [0,0];
                    post = [1,1];
                end
            end
        end
    end
end

function [] = register_position(position, input_dir, mip_input_dir, output_dir)
    config = get_config();
    positions = config('positions');
    cycles = config('cycles');
    types = config('types');

    mip = read_mips(mip_input_dir, position);
    ims = read_position(input_dir, position);

    [imy,imx] = size(mip{1});
    [refsize,prev,pos,post] = position_in_combined(position, imx, imy);

    [topos,frompos] = trans_to_position(pos, imx, imy);
    if ~isnan(prev)
        mip_prev = read_mips(mip_input_dir, position-1);
        [topos_prev,~] = trans_to_position(prev, imx, imy);
    end
    if ~isnan(post)
        mip_post = read_mips(mip_input_dir, position+1);
        [topos_post,~] = trans_to_position(post, imx, imy);
        
    end

    % now register the images with each other
    combined = cell([numel(cycles) 1]);
    for cycleidx=1:numel(cycles)
        cycle = cycles(cycleidx);
        im = mip{cycleidx};
        expanded = imwarp(im, affine2d(topos), 'cubic', 'OutputView', refsize);
        if config('debug') ; figure,imshow(expanded); title(['Expanded ' num2str(cycle)]); disp('Expanded:'); disp(topos); end;
        if ~isnan(prev)
            expanded_prev = imwarp(mip_prev{cycleidx}, affine2d(topos_prev), 'cubic', 'OutputView', refsize);
            if config('debug') ; figure,imshow(expanded_prev); title('Expanded prev'); disp('Expanded prev:'); disp(topos_prev); end;
            [slice,slice_prev,toslice_prev] = slice_position(pos, prev, expanded, expanded_prev, imx, imy);
            if config('debug_slices') ; figure,imshow(imfuse(slice,slice_prev,'montage')); title('Slice prev'); end
            % now register the slices
            %toreg_prev = register_images(slice_prev, slice);
            toreg_prev = register_with_gradient(slice_prev, slice);
            if config('debug') ; title(['Cycle ' num2str(cycle)]); end
            reg_prev = imwarp(mip_prev{cycleidx}, affine2d(toslice_prev * toreg_prev * itrans(toslice_prev) * topos_prev), 'cubic', 'OutputView', refsize);
        end
        if ~isnan(post)
            expanded_post = imwarp(mip_post{cycleidx}, affine2d(topos_post), 'cubic', 'OutputView', refsize);
            if config('debug') ; figure,imshow(expanded_post); title('Expanded post'); disp('Expanded post:'); disp(topos_post); end;
            [slice,slice_post,toslice_post] = slice_position(pos, post, expanded, expanded_post, imx, imy);
            if config('debug_slices') ; figure,imshow(imfuse(slice,slice_post,'montage')); title('Slice post');end
            % now register the slices
            %toreg_post = register_images(slice_post, slice);
            toreg_post = register_with_gradient(slice_post, slice);
            if config('debug') ; title(['Cycle ' num2str(cycle)]); end
            reg_post = imwarp(mip_post{cycleidx}, affine2d(toslice_post * toreg_post * itrans(toslice_post) * topos_post), 'cubic', 'OutputView', refsize);
        end

        % produce the combined image
        combined{cycleidx} = expanded;
        if ~isnan(prev)
            combined{cycleidx} = max(combined{cycleidx}, reg_prev);
        end
        if ~isnan(post)
            combined{cycleidx} = max(combined{cycleidx}, reg_post);
        end

        toreg = [1 0 0; 0 1 0; 0 0 1];
        if cycleidx > 1
            % register against the first combined image
            %toreg = register_images(combined{cycleidx}, combined{1});
            toreg = register_with_gradient(combined{cycleidx}, combined{1});
        end

        % transform all images for this cycle and the mip, then save them
        warped = imwarp(im, affine2d(topos * toreg * frompos), 'cubic', 'OutputView', imref2d(size(im)));
        save_mip(output_dir, position, cycle, warped);

        % save the individual images
        %registered = imwarp(ims{cycleidx}, affine2d(topos * toreg * frompos), 'cubic', 'OutputView', imref2d(size(im)));
        %save_images(output_dir, position, cycle, registered);
    end
end

function [] = save_images(output_dir, position, cycle, registered)
    config = get_config();
    cycles = config('cycles');
    types = config('types');

    % write out all images
    fmt = '%s/%s';

    % write the registered images
    for channel = 1:size(registered, 3)
        fn = sprintf(fmt, output_dir, sprintf(config('img_cycle'), position, cycle, types{channel}));
        write_image(registered(:,:,channel), fn);
    end
end

function [] = register_prev(position, reg_with_position, input_dir, mip_input_dir, output_dir)
    config = get_config();
    cycles = config('cycles');
    positions = config('positions');
    I = [1 0 0; 0 1 0; 0 0 1];
    
    if isempty(find(positions == reg_with_position,1))
        error(['Position ' num2str(reg_with_position) ' is not included in this data set.']);
    end
    
    reg = read_mips(mip_input_dir, position);
    regwith = read_mips(mip_input_dir, reg_with_position);
    
    [imy,imx] = size(reg{1});
    [toprev,fromprev, ...
     tocomb,fromcomb, ...
     refsize, ...
     slice_x1,slice_x2, ...
     slice_y1,slice_y2, ...
     slice_reg_x1,slice_reg_x2, ...
     slice_reg_y1,slice_reg_y2] = relative_to_prime(position, reg_with_position, imx, imy);
    
    % estimate the transform needed to register each pair
    % of images from the same cycle.
    before = cell([numel(cycles) 1]);
    after = cell([numel(cycles) 1]);
    for cycleidx=1:numel(cycles)
        cycle = cycles(cycleidx);
        disp(['Cycle ' num2str(cycle) '...']);
        
        im = reg{cycleidx};
        imr = regwith{cycleidx};
        
        slice = im(slice_y1:slice_y2,slice_x1:slice_x2);
        slice_reg = imr(slice_reg_y1:slice_reg_y2,slice_reg_x1:slice_reg_x2);
        toreg = register_images(slice, slice_reg);
        if config('debug') ; title(['Cycle ' num2str(cycle)]); end
        
        before{cycleidx} = I * toreg * toprev * tocomb;
        after{cycleidx} = fromcomb * fromprev * I;
    end
    
    % produce a combined image from the pairs in the same cycle
    for cycleidx=1:numel(cycles)
        im = reg{cycleidx};
        imr = regwith{cycleidx};
        cur = combine(im, imr, before{cycleidx}, tocomb, refsize);
        
        if cycleidx > 1
            % register against first combined image
            toreg = register_images(cur, first_combined);
            if config('debug') ; title(['Cycle ' num2str(cycle)]); end
            before{cycleidx} = before{cycleidx} * toreg;
        else 
            first_combined = cur;
        end
    end
    
    % do the full transformation
    for cycleidx=1:numel(cycles)
        cycle = cycles(cycleidx);
        % warp the mips
        im = reg{cycleidx};
        warped = imwarp(im, affine2d(before{cycleidx} * after{cycleidx}), 'cubic', 'OutputView', imref2d(size(im)));
        save_mip(output_dir, position, cycle, warped);
        % save each individual image
        % FIXME
        
    end
end

%{
function [] = register_cycles(input_dir, mip_input_dir, output_dir)
    config = get_config();
    write_config(output_dir);
    positions = config('positions');
    for posidx = 1:numel(positions)
        disp(['Position ' num2str(posidx) '...']);
        position = positions(posidx);
        
        if position == 1
            reg_with_position = 2;
        else
            reg_with_position = position - 1;
        end
        
        register_with(position, reg_with_position, input_dir, mip_input_dir, output_dir);
    end
end

function [] = foobar(input_dir, mip_input_dir, output_dir) 
% register_cycles
% register pairs of images, then register the pairs against pairs from 
% other cycles 

    config = get_config();
    write_config(output_dir);
    cycles = config('cycles');
    positions = config('positions');
    regby = config('register_by');
    cpd_opt = config('cpd_opt');
    
    % identity transforms for every cycle and position
    before = init_transforms();
    after = init_transforms();
    
    for posidx = 1:numel(positions)
        position = positions(posidx);
        disp(['Position ' num2str(position) '...']);
        
        % read the mips for this position, all cycles.
        m = read_mips(mip_input_dir, position);
        % get the position info for this position
        [rownum, rowpos, ~, first_in_row, ~, numcols] = position_info(posidx);
        
        % run registration for images which are not the first on a row
        % otherwise just shift the image
        run_registration = (strcmp(regby, 'row') && ~first_in_row) ...
                           || (~strcmp(regby, 'row') && posidx > 1);
        
        for cycleidx=1:numel(cycles)
            cycle = cycles(cycleidx);
            disp(['Cycle ' num2str(cycle) '...']);
            
            if ~run_registration
                % just shift the image, if necessary
                [~, ~, toprev, fromprev] = shift(posidx, m{cycle});
                before{cycleidx,posidx} = before{cycleidx,posidx} * toprev;
                after{cycleidx,posidx} = after{cycleidx,posidx} * fromprev;
                
                previous_unregistered = true;
                previous = m;
                continue;
            end
            
            % register each cycle image with the same cycle from the
            % previous position
            [toprev, toreg, fromprev] = register_with_previous(posidx, previous{cycle}, m{cycle});
            before{cycleidx,posidx} = before{cycleidx,posidx} * toreg * toprev;
            after{cycleidx,posidx} = after{cycleidx,posidx} * fromprev;
            
            
            % produce a combined pair of images to register against other
            % cycles for the same positions
            [vsize,hsize] = pair_bbox(m{cycle}) ;
            warped = imwarp(m{cycle}, affine2d(before{cycleidx,posidx}), ...
                'cubic', 'OutputView', imref2d([vsize,hsize]));
            fused = imfuse(previous{cycle}, warped, imref2d([vsize,hsize]));
            
            if config('debug')
                figure;
                imshow(imfuse(previous{cycle}, warped, 'ColorChannels', 'red-cyan'));
                title(['pos:' num2str(posidx)]);
            end
            
            % save the current position for the next round of registration
            previous = m;
        end
    end
end
%}

%{
function [] = register_with(position, reg_with_position, input_dir, mip_input_dir, output_dir)
    config = get_config();
    cycles = config('cycles');
    positions = config('positions');
    I = [1 0 0; 0 1 0; 0 0 1];
    
    if isempty(find(positions == reg_with_position,1))
        error(['Position ' num2str(reg_with_position) ' is not included in this data set.']);
    end
    
    reg = read_mips(mip_input_dir, position);
    regwith = read_mips(mip_input_dir, reg_with_position);
    
    [imy,imx] = size(reg{1});
    [toprev,fromprev, ...
     tocomb,fromcomb, ...
     refsize, ...
     slice_x1,slice_x2, ...
     slice_y1,slice_y2, ...
     slice_reg_x1,slice_reg_x2, ...
     slice_reg_y1,slice_reg_y2] = relative_to(position, reg_with_position, imx, imy);
    
    % estimate the transform needed to register each pair
    % of images from the same cycle.
    before = cell([numel(cycles) 1]);
    after = cell([numel(cycles) 1]);
    for cycleidx=1:numel(cycles)
        cycle = cycles(cycleidx);
        disp(['Cycle ' num2str(cycle) '...']);
        
        im = reg{cycleidx};
        imr = regwith{cycleidx};
        
        slice = im(slice_y1:slice_y2,slice_x1:slice_x2);
        slice_reg = imr(slice_reg_y1:slice_reg_y2,slice_reg_x1:slice_reg_x2);
        toreg = register_images(slice, slice_reg);
        if config('debug') ; title(['Cycle ' num2str(cycle)]); end
        
        before{cycleidx} = I * toreg * toprev * tocomb;
        after{cycleidx} = fromcomb * fromprev * I;
    end
    
    % produce a combined image from the pairs in the same cycle
    for cycleidx=1:numel(cycles)
        im = reg{cycleidx};
        imr = regwith{cycleidx};
        cur = combine(im, imr, before{cycleidx}, tocomb, refsize);
        
        if cycleidx > 1
            % register against first combined image
            toreg = register_images(cur, first_combined);
            if config('debug') ; title(['Cycle ' num2str(cycle)]); end
            before{cycleidx} = before{cycleidx} * toreg;
        else 
            first_combined = cur;
        end
    end
    
    % do the full transformation
    for cycleidx=1:numel(cycles)
        cycle = cycles(cycleidx);
        % warp the mips
        im = reg{cycleidx};
        warped = imwarp(im, affine2d(before{cycleidx} * after{cycleidx}), 'cubic', 'OutputView', imref2d(size(im)));
        save_mip(output_dir, position, cycle, warped);
        % save each individual image
        % FIXME
        
    end
end
%}

function [] = save_mip(output_dir, position, cycle, im)
    config = get_config();
    fmt = '%s/%s';
    fn = sprintf(fmt, output_dir, sprintf(config('img_mip'), position, cycle));
    write_image(im, fn);
end

function [combined] = combine(im, im_reg, before, tocomb, refsize)
    warped = imwarp(im, affine2d(before), 'cubic', 'OutputView', refsize);
    warped_reg = imwarp(im_reg, affine2d(tocomb), 'cubic', 'OutputView', refsize);
    % make, essentially, a MIP of these two.
    combined = max(warped,warped_reg);
end

function [toreg] = register_images(im, im_reg)
    config = get_config();
    thresh = config('surf_features_metric_threshold');
    points = detectSURFFeatures(im_reg, 'MetricThreshold', thresh);
    [features, valid] = extractFeatures(im_reg,  points);
    r_points = detectSURFFeatures(im, 'MetricThreshold', thresh);
    [r_features, r_valid] = extractFeatures(im, r_points);
    
    pairs = matchFeatures(features, r_features);
    matched  = valid(pairs(:,1));
    r_matched = r_valid(pairs(:,2));
    [tform,~,~] = estimateGeometricTransform(r_matched,matched,'affine');
    toreg = tform.T;
    
    if config('debug')
        szz = imref2d(size(im_reg));
        figure,imshow(...
            imfuse(im_reg, ...
                   imwarp(im, affine2d(toreg), 'cubic', 'OutputView', szz), ...
                   'ColorChannels', 'red-cyan'));
    end
    
    if config('registration_ignore_rotation')
        toreg = [1 0 0; 0 1 0; toreg(3,1) toreg(3,2) 1];
    end
end

function [toreg] = register_with_gradient(im, im_reg)
    config = get_config();
    % images produced on same device
    [optimizer, metric] = imregconfig('monomodal');
    optimizer.MaximumIterations = config('slice_registration_iterations');
    optimizer.MaximumStepLength = config('slice_registration_steplength_max');
    optimizer.MinimumStepLength = config('slice_registration_steplength_min');
    optimizer.RelaxationFactor = config('slice_registration_relaxation_factor');
    disp(['Registering images with ' num2str(optimizer.MaximumIterations) ' iterations']);
    tform = imregtform(im, im_reg, 'rigid', optimizer, metric);
    toreg = tform.T;
end

%{
function [regtoprev,fromprev, ...
          tocomb,fromcomb, ...
          refsize, ...
          slice_x1,slice_x2, ...
          slice_y1,slice_y2, ...
          slice_reg_x1,slice_reg_x2, ...
          slice_reg_y1,slice_reg_y2] = relative_to_prime(position, reg_with_position, imx, imy) 
end
%}

function [toprev,fromprev, ...
          tocomb,fromcomb, ...
          refsize, ...
          slice_x1,slice_x2, ...
          slice_y1,slice_y2, ...
          slice_reg_x1,slice_reg_x2, ...
          slice_reg_y1,slice_reg_y2] = relative_to(position, reg_with_position, imx, imy) 
    % FIXME this is not at all generic, it relies on the knowledge that
    % there are 16 positions and a zig-zag pattern
    config = get_config();
    if ~(reg_with_position == (position+1) || reg_with_position == (position-1))
        error('Only registration with the previous or next image is supported.');
    end
        
    [~, ~, oddrow, first_in_row, ~, ~] = position_info(position);
    [~, ~, ~, last_in_row, ~, ~] = position_info(reg_with_position);
    reg_with_previous = reg_with_position == (position-1);
    
    tocomb = [1 0 0; 0 1 0; 0 0 1];
    shifty = imy - config('position_overlap_v');
    shiftx = imx - config('position_overlap_h');
    vref = imref2d([((2*imy) - config('position_overlap_v')) imx]);
    href = imref2d([imy ((2*imx) - config('position_overlap_h'))]);
    if first_in_row && reg_with_previous
        % this image is below the image with which we will register
        toprev = [1 0 0; 0 1 0; 0 shifty 1];
        
        slice_x1 = 1; slice_x2 = imx;
        slice_y1 = 1; slice_y2 = config('position_overlap_v');
        slice_reg_x1 = 1; slice_reg_x2 = imx;
        slice_reg_y1 = imy - config('position_overlap_v'); slice_reg_y2 = imy;
        
        refsize = vref;
    elseif last_in_row && ~reg_with_previous
        % this image is above the image with which we will register
        toprev = [1 0 0; 0 1 0; 0 -shifty 1];
        tocomb = [1 0 0; 0 1 0; 0 shifty 1];
        
        slice_x1 = 1; slice_x2 = imx;
        slice_y1 = imy - config('position_overlap_v'); slice_y2 = imy;
        slice_reg_x1 = 1; slice_reg_x2 = imx;
        slice_reg_y1 = 1; slice_reg_y2 = config('position_overlap_v');
        
        refsize = vref;
    elseif (oddrow && reg_with_previous) || (~oddrow && ~reg_with_previous)
        % this image is to the left of the image with which we will register
        toprev = [1 0 0; 0 1 0; 0 -shiftx 1];
        tocomb = [1 0 0; 0 1 0; 0 shiftx 1];
        
        slice_x1 = imx - config('position_overlap_h'); slice_x2 = imx;
        slice_y1 = 1; slice_y2 = imy;
        slice_reg_x1 = 1; slice_reg_x2 = config('position_overlap_h');
        slice_reg_y1 = 1; slice_reg_y2 = imy;
        
        refsize = href;
    elseif (oddrow && ~reg_with_previous) || (~oddrow && reg_with_previous)
        % this image is to the right of the image with which we will register
        toprev = [1 0 0; 0 1 0; 0 shiftx 1];
        
        slice_x1 = 1; slice_x2 = config('position_overlap_h');
        slice_y1 = 1; slice_y2 = imy;
        slice_reg_x1 = imx - config('position_overlap_h'); slice_reg_x2 = imx;
        slice_reg_y1 = 1; slice_reg_y2 = imy;
        
        refsize = href;
    else
        error('Kellen, you missed something.');
    end
    fromprev = [1 0 0; 0 1 0; -toprev(3,1) -toprev(3,2) 1];
    fromcomb = [1 0 0; 0 1 0; -tocomb(3,1) -tocomb(3,2) 1];
end


function [] = foo(input_dir, mip_input_dir, output_dir) 
    config = get_config();
    write_config(output_dir);
    cycles = config('cycles');
    positions = config('positions');
    regby = config('register_by');
    cpd_opt = config('cpd_opt');

    points = cell([numel(cycles) 1]);
    ims = cell([numel(cycles) 1]);
    if strcmp(regby, 'row')
        [~, total_rows] = dims();
        points = cell([numel(cycles) total_rows]);
        ims = cell([numel(cycles) total_rows]);
    end
    
    before = init_transforms();
    after = init_transforms();
    
    for cycleidx=1:numel(cycles)
        cycle = cycles(cycleidx);
        disp(['Cycle ' num2str(cycle) '...']);
        m = read_cycle_mips(mip_input_dir, cycle);
        
        disp('Registering mips...');
        for posidx = 1:numel(positions) 
            position = positions(posidx);
            disp(['Position ' num2str(position) '...']);
            
            [rownum, rowpos, ~, first_in_row, ~, numcols] = position_info(posidx);
            
            % true = run registration, otherwise just shift the image
            normal_registration_pos = (strcmp(regby, 'row') && ~first_in_row) || (~strcmp(regby, 'row') && posidx > 1);
            
            % register/shift the image
            if normal_registration_pos
                [toprev, toreg, fromprev] = register_with_previous(posidx, previous,  m{position});
                before{cycleidx,posidx} = before{cycleidx,posidx} * toreg * toprev;
                after{cycleidx,posidx} = after{cycleidx,posidx} * fromprev;
            else
                [~, ~, toprev, fromprev] = shift(posidx, m{position});
                before{cycleidx,posidx} = before{cycleidx,posidx} * toprev;
                after{cycleidx,posidx} = after{cycleidx,posidx} * fromprev;
            end
            
            % produce the image for the next image to register against
            [vsize,hsize] = bbox(m{position}) ;
            warped = imwarp(m{position}, affine2d(before{cycleidx,posidx}), ...
                'cubic', 'OutputView', imref2d([vsize,hsize])); 
            
            if normal_registration_pos
                figure,imshow(imfuse(previous,warped, 'ColorChannels', 'red-cyan'));title(['pos:' num2str(posidx)]);
                % make a sort of MIP for these two images to fuse them
                previous = max(previous, warped);
            else
                previous = warped;
            end
            ims{cycleidx,rownum+1} = previous;
            
            if strcmp(regby, 'row') && numcols == (rowpos+1)
                points{cycleidx,rownum+1} = pointdrift_points(previous);
                if cycleidx > 1
                    % register against the first cycle
                    origpoints = points{1,rownum+1};
                    [trans,~] = cpd_register(origpoints, points{cycleidx,rownum+1}, cpd_opt);
                    tform = affine2d([trans.R [0;0]; [trans.t' 1]]);
                    
                    firstim = ims{1,rownum+1};
                    figure,imshow(imfuse(firstim, ...
                                         imwarp(ims{cycleidx,rownum+1}, ...
                                                 tform, ...
                                                 'cubic', 'OutputView', ...
                                                 imref2d(size(firstim)))));
                    title(['cycle:' num2str(cycleidx)])
                end
            end
            
            %figure,imshow(previous);title(['pos:' num2str(posidx)]);
        end
        figure,imshow(previous);
        keyboard;

        %{
        disp('Finding points...');
        % get the interesting points for the entire fucking thing
        points{cycleidx} = pointdrift_points(bigfuckingimage);
        %}
    end
    
    % FIXME now register the points and record the shit in the transforms
    % then apply to the individual images
end

function [rownum, rowpos, oddrow, first_in_row, numrows, numcols, last_in_row] = position_info(posidx) 
    % zig zag makes this weird; these numbers only work for tilings which
    % have an even number of tiles. there is likely a better way to do this
    % more generically...
    posidx = posidx - 1;
    
    [numcols, numrows] = dims();
    rownum = fix(posidx / numcols);
    rowpos = mod(posidx, numcols);
    first_in_row = 0 == rowpos;
    last_in_row = rowpos == (numcols-1);
    oddrow = (1 == mod(rownum, 2));
end

function [shiftx, shifty, toprev, fromprev] = shift(posidx, im)
    config = get_config();
    regby = config('register_by');
    
    [imy,imx] = size(im);
    [rownum, rowpos, oddrow, ~, ~, numcols] = position_info(posidx);
    
    shiftx = (imx * rowpos) - (rowpos * config('position_overlap_h'));
    if oddrow
        shiftx = (imx * (numcols-1)) - (rowpos * imx) - (config('position_overlap_h') * (numcols-rowpos-1));
    end
    
    shifty = (imy * rownum) - (config('position_overlap_v') * rownum);
    if strcmp(regby, 'row')
        shifty = 0;
    end
    toprev = [1 0 0; 0 1 0; shiftx shifty 1];
    fromprev = [1 0 0; 0 1 0; -toprev(3,1) -toprev(3,2) 1];
end

function [numcols, numrows] = dims()
    config = get_config();
    numcols = config('position_row_length');
    numrows = ceil(max(config('positions')) / numcols);
end

function [vsize,hsize] = pair_bbox(im) 
    % produce a bounding box for pair-based registration L-R only
    s = [1 2] .* size(im);
    vsize = s(1);
    hsize = s(2);
end

function [vsize,hsize] = bbox(im) 
    config = get_config();
    regby = config('register_by');
    
    [numcols, numrows] = dims();
    [imy,imx] = size(im);
    vsize = (imy * numrows);
    if strcmp(regby, 'row')
        vsize = imy;
    end
    hsize = (imx * numcols);
end

function [toworld, fromworld] = world_space(posidx, im)
    config = get_config();
    [rownum, rowpos, oddrow, first_in_row, numrows, numcols] = position_info(posidx);
    [vsize,hsize] = bbox(im);
    I = [1 0 0; 0 1 0; 0 0 1];
    [imy,imx] = size(im);

    toworld = I;
    fromworld = I;
    toworld(3,1) = (imx * rowpos) - (rowpos * config('position_overlap_h'));
    %FIXME this part has to be relative to the last positions...?
    if oddrow
        toworld(3,1) = (hsize - (imx * (rowpos+1))) - ((numcols-rowpos-1) * config('position_overlap_h'));
    end
    toworld(3,2) = (imy * (rownum - 1)) - ((rownum-1) * config('position_overlap_v'));
    fromworld(3,1) = -toworld(3,1);
    fromworld(3,2) = -toworld(3,2); 
end

function [transforms] = init_transforms() 
    config = get_config();
    cycles = config('cycles');
    positions = config('positions');
    transforms = cell([numel(cycles) numel(positions)]);
    
    I = [1 0 0; 0 1 0; 0 0 1];
    for cycleidx=1:numel(cycles)
        for posidx = 1:numel(positions) 
            transforms{cycleidx,posidx} = I;
        end
    end
end

function [toprev, toreg, fromprev] = register_with_previous(posidx, previous, current)
    config = get_config();
    [my, mx] = size(current);
    [~, ~, oddrow, first_in_row, ~, ~] = position_info(posidx);
    [shiftx, shifty, toprev, fromprev] = shift(posidx, current);

    % get a strip of the previous image against which we will register
    if first_in_row
        % register horizontal strip
        prevslice = previous(1+shifty:(shifty+config('position_overlap_v')),1+shiftx:(shiftx+mx));
        curslice = current(1:config('position_overlap_v'),:);
    else
        % register vertical strip
        if oddrow
            prevslice = previous(1+shifty:(shifty+my),(1+shiftx+mx-config('position_overlap_h')):(shiftx+mx));
            curslice = current(:,end-config('position_overlap_h'):end);
        else
            prevslice = previous(1+shifty:(shifty+my),1+shiftx:(shiftx+config('position_overlap_h')));
            curslice = current(:,1:config('position_overlap_h'));
        end
        sz = imref2d(size(previous));
        %figure,imshow(imfuse(previous,imwarp(current, affine2d(toprev),'cubic', 'OutputView', sz), 'ColorChannels', 'red-cyan'));
        %keyboard;
    end
    
    points = detectSURFFeatures(prevslice, 'MetricThreshold', config('surf_features_metric_threshold'));
    [features, valid] = extractFeatures(prevslice,  points);
    r_points = detectSURFFeatures(curslice, 'MetricThreshold', config('surf_features_metric_threshold'));
    [r_features, r_valid] = extractFeatures(curslice, r_points);
    
    pairs = matchFeatures(features, r_features);
    matched  = valid(pairs(:,1));
    r_matched = r_valid(pairs(:,2));
    [tform,~,~] = estimateGeometricTransform(r_matched,matched,'affine');
    
    toreg = tform.T;
    if config('registration_ignore_rotation')
        toreg = [1 0 0 ; 0 1 0; toreg(3,1) toreg(3,2) 1];
    end
    
    if config('debug')
        szz = imref2d([size(prevslice,1) size(prevslice,2)*2]);
        if first_in_row
            szz = imref2d([2*size(prevslice,1) size(prevslice,2)]);
        end
        figure,imshow(imfuse(prevslice, ...
        imwarp(curslice, affine2d(toreg), 'cubic', 'OutputView', szz), 'ColorChannels', 'red-cyan'));
        title(['reg ' num2str(posidx)]);
    end

    %figure,imshow(imfuse(previous,imwarp(current, affine2d(toreg * toprev),'cubic', 'OutputView', sz), 'ColorChannels', 'red-cyan'));title('reggy');
    %keyboard;
    
    % if this is an image registered along the right hand side
    if oddrow && ~first_in_row
        % first shift the image to the correct location, then apply the
        % registration transform, then shift back
        toreg = [1 0 0; 0 1 0; -(mx-config('position_overlap_h')) 0 1] * toreg * [1 0 0 ; 0 1 0 ; (mx-config('position_overlap_h')) 0 1];
    end
    disp('Reg:');disp(toreg);
    disp('To:');disp(toprev);
    disp('From:');disp(fromprev);
end
