function [r] = roi(im) 
    config = get_config();
    r = im(config('roi_y'),config('roi_x'),:);
end