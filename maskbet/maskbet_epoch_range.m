function epochRange = maskbet_epoch_range(epochName, plotOrAnalyze)

% Returns a range of time around the alignement time for an input epoch


% Default is to return an analysis epoch
if nargin < 2
    plotOrAnalyze = 'analyze';
end



switch plotOrAnalyze
    case 'plot'
        % Set a standard range so each panel of a graph
        epochRange = -349 : 350;
        switch epochName
            case 'decFixWindowEntered'
                shiftRange = 100;
            case 'decTargOn'
                shiftRange = 100;
            case 'decResponseOnset'
                shiftRange = 0;
            case 'betFixWindowEntered'
                shiftRange = 100;
            case 'betTargOn'
                shiftRange = 100;
            case 'betResponseOnset'
                shiftRange = 0;
            case 'rewardOn'
                shiftRange = 100;
        end
        
        
        
        epochRange = epochRange + shiftRange;
    case 'analyze'
end
