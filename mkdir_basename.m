function [ ] = mkdir_basename(fn)
    [dir,~,~] = fileparts(fn);
    [s,mess,~] = mkdir(dir);
    if ~s
        error(['Could not create ' dir ' Message:' mess]);
    end
end