function [ L ] = wshed( bw )
%WSHED perform watershed on distance transform of bw

    D = bwdist(~bw);
    D = -D;
    D(~bw) = -Inf;
    L = watershed(D);
    % find the background label and remove the bg class
    [r, c] = find(bw == 0, 1);
    bglabel = L(r, c);
    L(L == bglabel) = 0;
end