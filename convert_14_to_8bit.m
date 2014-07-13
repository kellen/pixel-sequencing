function [ ] = convert_14_to_8bit( input_path, output_path )
%CONVERT_14_TO_8BIT Finds all folders 1 level under input_path, descends
%into them and converts the images to 8-bit images on the assumption that 
%they're 14-bit images in 16-bit containers and puts the results in 
%output_path. Terrible.

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
                % convert it, assuming 16-bit file  and 14-bit actual values
                i2 = uint8(floor(255 * (double(i) / (2^14))));
                
                % save it
                imwrite(i2,[output_path '/' d.name '/' f.name]);
            end      
        end
    end
end

