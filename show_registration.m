function [] = show_registration(before_input_dir, after_input_dir, before_do_input_dir, position, output_dir) 

    before = read_mips(before_input_dir, position);
    make_color_roi_comparision(before, sprintf('%s/%s', output_dir, 'reg-before.tif'));
    after = read_mips(after_input_dir, position);
    make_color_roi_comparision(after, sprintf('%s/%s', output_dir, 'reg-after.tif'));
    
    total = max(after{1}, after{2});
    total = max(total, after{3});
    total = max(total, after{4});
    
    pos = read_position(before_do_input_dir, position);
    before_do = pos('do');
    pos = read_position(after_input_dir, position);
    after_do = pos('do');
    
    write_image(roi(imfuse(total, before_do)), sprintf('%s/%s', output_dir, 'reg-do-before.tif'));
    write_image(roi(imfuse(total, after_do)), sprintf('%s/%s', output_dir, 'reg-do-after.tif'));
end

function [] = make_color_roi_comparision(mips, fn)
    first = imfuse(mips{1}, mips{2}, 'ColorChannels', 'red-cyan');
    %figure, imshow(first);
    second = imfuse(mips{3}, mips{4}, 'ColorChannels', 'green-magenta');
    %figure, imshow(second);
    final = imfuse(first, second, 'blend');
    %figure, imshow(final);
    
    write_image(roi(final), fn);
end