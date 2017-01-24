function epochRange = ccm_epoch_range(epochName, breadth, rangeFactor)

% Returns a range of time around the alignement time for an input epoch

% plotOrAnalyze: input: 'plot', 'analyze', 'broad'.  Might want different window depending on whether we are
% viewing the data or analyzing a specific epoch.

% rangeFactor: 'true' or 'false'. Doubles the normal range. Might want to
% see out further than hard-coded here.

% Default is to return an analysis epoch
if nargin < 2
    breadth = 'analyze';
end
if nargin < 3
    rangeFactor = 1;
end



switch breadth
    case 'plot'
        % Set a standard range so each panel of a graph
        epochRange = -99 : 400;
        epochRange = -99 : 300;
        switch epochName
            case 'fixWindowEntered'
                shiftRange = 0;
            case 'targOn'
                shiftRange = 0;
            case 'checkerOn'
                shiftRange = 50;
            case 'stopSignalOn'
                shiftRange = -100;
            case 'responseOnset'
                shiftRange = -100;
            case 'responseOffset'
                shiftRange = -150;
            case 'toneOn'
                shiftRange = -50;
            case 'rewardOn'
                shiftRange = -50;
        end
        
        
        
    case 'analyze'
        % Set a broad range for analyses
        epochRange = -299 : 500;
                shiftRange = 0;
%         switch epochName
%             case 'fixWindowEntered'
%                 shiftRange = -100;
%             case 'targOn'
%                 shiftRange = 0;
%             case 'checkerOn'
%                 shiftRange = 50;
%             case 'stopSignalOn'
%                 shiftRange = -100;
%             case 'responseOnset'
%                 shiftRange = -100;
%             case 'responseOffset'
%                 shiftRange = -150;
%             case 'toneOn'
%                 shiftRange = -50;
%             case 'rewardOn'
%                 shiftRange = -50;
%         end
    case 'broad'
        epochRange = -99 : 300;
end
        epochRange = epochRange + shiftRange;

epochRange = epochRange * rangeFactor;
