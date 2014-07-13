function [ ] = write_image(img, fn)
    config = get_config();
    filetype = config('filetype');
    mkdir_basename(fn);
    imwrite(img, fn, filetype);
end