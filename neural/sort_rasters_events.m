function [alignedrasters, event_marks] = sort_rasters_events(alignedrasters, event_marks)

% function [alignedrasters, event_marks] = sort_rasters_events(alignedrasters, event_marks)
%
% Sorts a set of rasters in order of the event_marks that accompany them,
% earliest event_marks on bottom row to latest event marks on top row
% (w.r.t. alignment point).
%
% pgm: 09-2011

if ~isempty(alignedrasters)
    if isempty(event_marks)
        event_marks = [];
    else
        to_sort = [event_marks, alignedrasters];
        sorted = sortrows(to_sort, 1);
        alignedrasters = sorted(:, 2 : end);
        event_marks = sorted(:, 1);
    end
else
    alignedrasters = [];
    event_marks = [];
    
end