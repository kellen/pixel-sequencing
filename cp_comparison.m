function [] = cp_comparison(image_input_dir, mip_input_dir, centroids_input_dir, cp_input_dir, output_dir) 
    config = get_config();
    write_config(output_dir);

    positions = config('positions'); 
    valid = config('valid');
    blob_quality_threshold = config('threshold_quality_cellprofiler');
    
    fmt = '%s/%s';
    for posidx = 1:numel(positions)
        position = positions(posidx);
        disp(['Position ' num2str(position) '...']);
        
        if ismember(position, config('skip'))
            disp('Skipping...');
            continue;
        end
        
        pos = read_position(image_input_dir, position);
        mips = read_mips(mip_input_dir, position);
        
        mmip = max(mips{1}, mips{2});
        mmip = max(mmip, mips{3});
        mmip = max(mmip, mips{4});
        
        % load blobs
        [blobpositions, bloblabels, blobquality] = cp_centroids(cp_input_dir, position);
        % filter the blob positions based on the config
        blobpositions = blobpositions(blobquality > blob_quality_threshold, :);
        bloblabels = bloblabels(blobquality > blob_quality_threshold);
        
        % load centroids
        fn = sprintf(fmt, centroids_input_dir, sprintf(config('mat_centroids'), position));
        S = load(fn);
        centpositions = S.centpositions;
        centlabels = S.centlabels;
        
        % make figures
        valid_cent = ismember(centlabels, valid);
        valid_blob = ismember(bloblabels, valid);
        cent_precision = sum(valid_cent) / size(centlabels, 1);
        blob_precision = sum(valid_blob) / size(bloblabels, 1);

        blob_ok = blobpositions(valid_blob, :);
        blob_err = blobpositions(~valid_blob, :);
        cent_ok = centpositions(valid_cent, :);
        cent_err = centpositions(~valid_cent, :);
        
        disp(num2str(size(centlabels, 1) - sum(valid_cent)));
        
        %disp(['CP precision: ' num2str(blob_precision) ' (' num2str(sum(valid_blob)) '/' num2str(size(bloblabels, 1)) ')' ...    
        %' Per-pixel precision: ' num2str(cent_precision) ' (' num2str(sum(valid_cent)) '/' num2str(size(centlabels, 1)) ')']);
       
        if position < 2
        % bad colors, but what are good ones?
        sz = 20;
        blue = [0,0,1];
        red = [1,0,0];
        teal = [0,1,1];
        mag = [1,0,1];
        figure;
        imshow(mmip, [0 255]);
        hold on, scatter(blob_ok(:, 1), blob_ok(:, 2), sz, repmat(blue,size(blob_ok,1),1), '+'),
        hold on, scatter(blob_err(:, 1), blob_err(:, 2), sz, repmat(mag,size(blob_err,1),1), '+'),
        hold on, scatter(cent_ok(:, 1), cent_ok(:, 2), sz, repmat(blue,size(cent_ok,1),1), 'o'),
        hold on, scatter(cent_err(:, 1), cent_err(:, 2), sz, repmat(mag,size(cent_err,1),1), 'o');
    
        % do it the slow way
        % should probably accumulate all points to be drawn
        % in separate matrixes then draw one time rather than for each pair
        [closest_indexes, closest_dists] = dsearchn(blobpositions, delaunayn(blobpositions), centpositions);
        un = unique(closest_indexes);
        [counts, ~] = histc(closest_indexes,un);

        ign_err = 0;
        ign_ok = 0;
        same = 0;
        diff = 0;
        avg_same = 0;
        avg_diff = 0;
        for i=1:size(closest_indexes, 1)
            hold on;
            to_idx = closest_indexes(i);
            to = blobpositions(to_idx, :);
            c_idx = find(un == to_idx);

            % if there is more than one cent pointing to the same blob, only 
            % draw a line for the closest of them
            if counts(c_idx) > 1 && closest_dists(i) ~= min(closest_dists(closest_indexes == to_idx))
                % in this case, re-draw the glyph with a different color
                if ismember(centlabels(i), valid)
                    c = 'g';
                    ign_ok = ign_ok + 1;
                else
                    c = 'r';
                    ign_err = ign_err + 1;
                end
                scatter(centpositions(i, 1), centpositions(i, 2), sz, c, 'o');
                % continue with out drawing a connecting line
                continue;
            end

            % color the lines according to the correspondence between source
            % and destination labels

            if bloblabels(to_idx) == centlabels(i)
                c = 'g';
                same = same + 1;
                avg_same = avg_same + closest_dists(i);
            else
                c = 'r';
                diff = diff + 1;
                avg_diff = avg_diff + closest_dists(i);
            end

            from = centpositions(i, :);
            plot([from(1,1) to(1,1)], [from(1,2) to(1,2)], c);
        end
        avg_total = (avg_diff + avg_same)/(diff + same);
        avg_diff = avg_diff/diff;
        avg_same = avg_same/same;
        
        keyboard;

 disp(['Counts: same: ' num2str(same) ' (' num2str(same/(same + diff)) ') different: ' num2str(diff) ' ignored ok: ' num2str(ign_ok) ' ignored err: ' num2str(ign_err)]);
        disp(['Average distances: same: ' num2str(avg_same) ' diff: ' num2str(avg_diff) ' total: ' num2str(avg_total)]);
        end
    end
end