function [types] = types_list() 
    % T=1, G=2, C=3, A=4 to correspond with ID_list_BCpanel
    % this is the order in which the image channels come, 
    % i.e. T is the first channel (channel 2), G is the 2nd (channel 3)
    types = {'T', 'G', 'C', 'A'}; 
end