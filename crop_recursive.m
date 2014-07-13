function [  ] = crop_recursive( input_path, output_path )
%crop_recursive Recursively descends in input_path and crops all images
% to the specified size, writes them to output_path in the same structure.

    % make the output path
    [s,mess,~] = mkdir(output_path);
    if ~s
        error(['Could not create ' output_path ' Message:' mess]);
    end

    listing = dir(input_path);
    for d = listing'
        if d.isdir
            % make the subdir
            [s,mess,~] = mkdir(output_path, d.name);
            if ~s
                error(['Could not create ' d.name ' in' output_path ' Message:' mess]);
            end
            % find the images
            image_files = dir([input_path '/' d.name '/*.tif']);
            for f = image_files'
                % load each image
                i = imread([input_path '/' d.name '/' f.name]);
                
                % crop it, based on
                % magic numbers from image inspection
                i2 = i(55:1075, 25:1400, :);
                
                % save it
                imwrite(i2,[output_path '/' d.name '/' f.name]);
            end      
        end
    end

end

