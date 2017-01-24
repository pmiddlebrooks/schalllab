function [goFastTrial, goSlowTrial] = latency_match(goRT, Opt)
%
% function latency_match(goRT, stopRT, Opt)
% Opt: options structure contains:
%
%       matchMethod: 'ssrt' (default): match RTs based on whether RT is greater or less than SSRT + SSD
%                   'match': use a nearest neighbor algorithm and deliver
%                   equal number of go and stop trials
%       ssrt:   stop signal reaction time
%       ssd: stop signal delay

if nargin < 2
    Opt.matchMethod = 'ssrt';
    Opt.stopRT = []; % list of stop RTs ('match' method only)
    Opt.ssrt = [];   %  stop signal reaction time
    Opt.ssd = [];   %  stop signal delay
end

goFastTrial = [];
goSlowTrial = [];
if isempty(goRT)
    return
end

        stopLatency = Opt.ssd + Opt.ssrt;

switch Opt.matchMethod
    case 'ssrt'
        goSlowTrial = goRT > stopLatency;
        goFastTrial = goRT < stopLatency;
    case 'match'
        % Use nearest neighbor method to get rt-matched
        % trials
        error('Need to implement match method. Right now it uses knnsearch, whisch will return repeated latency-matched trial numbers from the pool of goRTs')
        % match the goRT and stopRT
        goFastTrial = knnsearch(goRT, stopRT);
        goSlowTrial = goRT > stopLatency;
        
        
%         nStopCorrect = size(iStopStopChecker.signal, 1);
%         data = ccm_match_rt(iGoTargChecker.eventLatency, iStopTargChecker.eventLatency(iStopTargTrial), nStopCorrect);
%         iGoSlowTrial = data.goSlowTrial;
%         
%         
%         
%         [rt, ind] = sort(goRT);
%         
%         lastInd = 1;
%         while nanmean(rt(lastInd:end)) < nanmean(stopRT)
%             goFastTrial = [goFastTrial; ind(lastInd)];
%             lastInd = lastInd + 1;
%             
%             
%         end
%         goSlowTrial = ind(lastInd:end);
%         %         trialList = 1 : length(goRT);
%         %         remaining = 1 : length(goRT);
%         %         remove = [];
%         %         goFastTrial = remaining;
%         %         while nanmean(goRT(remaining)) > nanmean(stopRT)
%         %             [y,i] = max(goRT(remaining));
%         %             remove = [remove; i]
%         %             goSlowTrial = [goSlowTrial; i];
%         %
%         %
%         %         end
%         %         goFastTrial = remaining;
        
        
end
