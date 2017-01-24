% Returns a defective CDF
% Normalized not to 1 but to response probability (or, overall error rate)
%
% RPH
function [CDF] = getDefectiveCDF(trialCorr, trialErr, SRT, plotFlag, ntiles)

zTrans = 0;
%if zTrans; disp('Computing Z-scores'); end


if nargin < 5; ntiles = .1:.1:.9; end
if nargin < 4; plotFlag = 0; end
if nargin < 3
    SRT = evalin('caller','SRT');
end

if isempty(trialCorr) && isempty(trialErr)
    CDF.correct = [NaN NaN];
    CDF.err = [NaN NaN];
    return
elseif isempty(trialCorr)
    errSortRT = sort(SRT(trialErr,1));
    errSortRT(find(isnan(errSortRT))) = [];
    resp_prob_err = 1;
    indices_err = ceil((ntiles) .* length(errSortRT));
    bins_err = errSortRT(indices_err);
    y_axis_err = (ntiles) .* resp_prob_err;
    
    CDF.correct = [NaN NaN];
    CDF.err(:,1) = bins_err;
    CDF.err(:,2) = y_axis_err;
elseif isempty(trialErr)
    corrSortRT = sort(SRT(trialCorr,1));
    corrSortRT(find(isnan(corrSortRT))) = [];
    resp_prob_correct = 1;
    indices_correct = ceil((ntiles) .* length(corrSortRT));
    bins_correct = corrSortRT(indices_correct);
    y_axis_correct = (ntiles) .* resp_prob_correct;
    
    CDF.correct(:,1) = bins_correct;
    CDF.correct(:,2) = y_axis_correct;
    CDF.err = [NaN NaN];
else
    if zTrans
        corrSortRT = sort(SRT(trialCorr,1));
        errSortRT = sort(SRT(trialErr,1));
        
        corrSortRT(find(isnan(corrSortRT))) = [];
        errSortRT(find(isnan(errSortRT))) = [];
        
        corrSortRT = zscore(corrSortRT,1);
        errSortRT = zscore(errSortRT,1);
    else
        corrSortRT = sort(SRT(trialCorr,1));
        errSortRT = sort(SRT(trialErr,1));
        
        corrSortRT(find(isnan(corrSortRT))) = [];
        errSortRT(find(isnan(errSortRT))) = [];
    end
    resp_prob_correct = length(trialCorr) / ( length(trialCorr) + length(trialErr) );
    resp_prob_err = 1 - resp_prob_correct;
    
    
    indices_correct = ceil((ntiles) .* length(corrSortRT));
    indices_err = ceil((ntiles) .* length(errSortRT));
    
    % indices_correct = [1 ceil((.1:.1:1) .* length(corrSortRT))];
    % indices_err = [1 ceil((.1:.1:1) .* length(errSortRT))];
    
    bins_correct = corrSortRT(indices_correct);
    bins_err = errSortRT(indices_err);
    
    
    % Accuracy rates are normalized to response probability
    
    y_axis_correct = (ntiles) .* resp_prob_correct;
    y_axis_err = (ntiles) .* resp_prob_err;
    
    % y_axis_correct = [0 (.1:.1:1) .* resp_prob_correct];
    % y_axis_err = [0 (.1:.1:1) .* resp_prob_err];
    
    
    CDF.correct(:,1) = bins_correct;
    CDF.correct(:,2) = y_axis_correct;
    CDF.err(:,1) = bins_err;
    CDF.err(:,2) = y_axis_err;
end
if plotFlag
    figure
    plot(bins_correct,y_axis_correct,'-ok',bins_err,y_axis_err,'--ok')
    fon
    title('Defective CDF')
end