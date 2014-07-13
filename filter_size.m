function [ newseq ] = filter_size(seq, size_threshold )
% removes objects smaller than size_threshold
    labels = unique(seq);
    labels = labels(labels ~= 0);
    newseq = zeros(size(seq));
    for i=1:numel(labels)
        bw = seq == labels(i);
        bw = bwareaopen(bw, size_threshold);
        newseq(bw) = labels(i);
    end
end