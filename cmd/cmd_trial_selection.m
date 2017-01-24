function  trialList = cmd_trial_selection(trialData, options)
% function  trialList = cmd_trial_selection(trialData, options)


if nargin < 2
   options.outcomeArray = 'all';
   options.ssdRange = 'all';
   options.targAngle = 'all';
   options.targHemifield = 'all';
   if nargin == 0
      trialList = options;
      return
   end
end


nTrial = size(trialData.trialOutcome, 1);
trialLogical = ones(nTrial, 1);






% Get list(s) of trials w.r.t. the outcome
if isfield(options, 'outcomeArray')
   if strcmp(options.outcomeArray, 'all')
      outcomeList = ones(nTrial, 1);
   else
      outcomeList = zeros(nTrial, 1);
      for iOutcomeIndex = 1 : length(options.outcomeArray)
         iOutcome = options.outcomeArray{iOutcomeIndex};
         
         outcomeList = outcomeList + strcmp(trialData.trialOutcome, iOutcome);
      end
   end
   trialLogical = trialLogical & outcomeList;
end



% Get list(s) of trials w.r.t. the SSDs
if isfield(options, 'ssdRange')
   if strcmp(options.ssdRange, 'none')
      % take any trials without a stop signal (nan values for ssd)
      trialLogical = trialLogical & isnan(trialData.ssd);
   elseif strcmp(options.ssdRange, 'any')
      % Do nothing- might want trials without regard to stop/go
   elseif strcmp(options.ssdRange, 'all')
      ssd = min(trialData.ssd) : max(trialData.ssd);
      trialLogical = trialLogical & ismember(trialData.ssd, ssd);
   else
      %     ssd = [options.ssdRange - 1, options.ssdRange, options.ssdRange + 1];
      % For now, use range to within an extra frame refresh
      ssd = [options.ssdRange - 13 : options.ssdRange + 18];
      trialLogical = trialLogical & ismember(trialData.ssd, ssd);
   end
end




% w.r.t. target hemifield
if isfield(options, 'targHemifield')
   targAngle = trialData.targAngle;
   if strcmp(options.targHemifield, 'all')
      % do nothing
   elseif strcmp(options.targHemifield, 'right')
      targTrial = ((targAngle > 270) & (targAngle <= 360)) | ...
         ((targAngle >= 0) & (targAngle <= 90)) | ...
         ((targAngle > -90) & (targAngle <= 0)) | ...
         ((targAngle >= -360) & (targAngle <= -270));
      trialLogical = trialLogical & targTrial;
   elseif strcmp(options.targHemifield, 'left')
      targTrial = ((targAngle > 90) & (targAngle <= 270)) | ...
         ((targAngle <= -90) & (targAngle > -270));
      trialLogical = trialLogical & targTrial;
   end
end



% w.r.t target angle
if isfield(options, 'targAngle')
   if strcmp(options.targAngle, 'all')
      % do nothing
   else
      trialLogical = trialLogical & ismember(trialData.targAngle, options.targAngle);
   end
end


trialList = find(trialLogical);
