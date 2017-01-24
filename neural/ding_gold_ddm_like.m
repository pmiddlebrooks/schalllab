function [ddmData] = ding_gold_ddm_like(data, varargin)

%

% Default analysis is done with neuronal spikes
dataType = 'spikes';
for i = 1 : 2 : length(varargin)
   switch varargin{i}
      case 'dataType'
         dataType = lower(varargin{i+1});
      otherwise
   end
end


% Unpac the data struct
switch dataType
   case 'spikes'
      trialMetric       = data.spikeRate;
      alignedRasters  = data.alignedRasters;
   case 'eeg'
      %       % Remove any data points with NaN for the trial Metric (happens when
      %       % inputing an eeg signal, e.g., with rejected trials (data switched to NaN)
      %       data(isnan(data.eegMeanEpoch), :) = [];
      trialMetric       = data.eegMeanEpoch;
      alignedSignal = data.alignedSignal;
   otherwise
      fprintf('ding_gold_ddm_like: invalid dataType entry: ''%s''/n', dataType);
      return
end
leftTrial       = data.leftTrial;
rightTrial      = data.rightTrial;
signalP         = data.signalP;
alignIndex      = data.alignIndex;

epochOffset     = data.epochOffset;


nTrial          = length(trialMetric);


% ________________________________________________________________
% DETERMINE WHETHER THE NEURON WAS "DDM-LIKE" ACCORDING TO DING & GOLD
% CRITERIA

% Initialize choice and coherence dependence to null assumption (i.e. they
% aren't)
choiceDependent     = false;
coherenceDependent  = false;
ddmLike             = false;

% some constants
alphaChoice     = .05;   % alpha criteria for choice dependence
alphaCoherence  = .05;   % alpha criteria for coherence dependence










% Choice dependence
leftMetric   = trialMetric(leftTrial);
rightMetric  = trialMetric(rightTrial);

[p, h, stats]   = ranksum(leftMetric, rightMetric);

if p < alphaChoice
   choiceDependent = true;
end

if nanmedian(leftMetric) >= nanmedian(rightMetric)
   leftIsIn    = true;
   inTrial     = leftTrial;
   outTrial    = rightTrial;
else
   leftIsIn    = false;
   inTrial     = rightTrial;
   outTrial    = leftTrial;
end








% Coherence dependence

% For IN trials


% Regress spikeRate vs signalStrength into RF
[coeffIn, sIn]          = polyfit(signalP(inTrial), trialMetric(inTrial), 1);
[yPredIn, deltaIn]  = polyval(coeffIn, signalP(inTrial), sIn);
statsIn             = regstats(signalP(inTrial), trialMetric(inTrial));
rIn                 = corr(signalP(inTrial), trialMetric(inTrial));

slopeIn     = coeffIn(1);
signSlopeIn = sign(slopeIn);
fTestIn     = statsIn.fstat.f;
pValIn      = statsIn.fstat.pval;


% Regress spikeRate vs signalStrength out of RF
[coeffOut, sOut]        = polyfit(signalP(outTrial), trialMetric(outTrial), 1);
[yPredOut, deltaOut] = polyval(coeffOut, signalP(outTrial), sOut);
statsOut            = regstats(signalP(outTrial), trialMetric(outTrial));
rOut                = corr(signalP(outTrial), trialMetric(outTrial));

slopeOut    = coeffOut(1);
signSlopeOut = -sign(slopeOut);  % NOTE THIS IS NEGATED BECAUSE IF IN CONDITION GOES HARD TO EASY, OUT GOES EASY TO HARD, VICE VERSA
fTestOut    = statsOut.fstat.f;
pValOut     = statsOut.fstat.pval;


% Decision tree to determine whether the neuron/signal was "coherence dependent"
if pValIn < alphaCoherence
   if pValOut > alphaCoherence
      coherenceDependent = true;
   elseif pValOut < alphaCoherence
      % slopeOut must have opposite sign than slopeIn
      if signSlopeIn ~= signSlopeOut
         coherenceDependent = true;
      end
   end
elseif pValIn > alphaCoherence
   if pValOut < alphaCoherence
      coherenceDependent = true;
   elseif pValOut > alphaCoherence
      % slopeOut must have opposite sign than slopeIn
      if signSlopeIn ~= signSlopeOut
         coherenceDependent = true;
      end
   end
end

if choiceDependent && coherenceDependent
   ddmLike = true;
end
% choiceDependent
% coherenceDependent
% ddmLike
%
%







% ________________________________________________________________
% DETERMINE THE ONSET OF THE CHOICE DEPENDENCE
tChoice = nan;  % Initialize to nan; in the case of eeg signals there may not be a tChoice (noise in signal)
nanThreshold = .3; % Establish a threshold of non-signal fraction over which we cease looking for tChoice
if ddmLike
   slideWindowWidth = 50; % ms, Ding and Gold used a 100ms sliding window
   slideWindowStep = 10; % ms, D&G used 10 ms steps
   
   
   choiceDependenceFound = false;
   iStepInd = 0;
   while ~choiceDependenceFound
      
      
      iEpochBegin = (alignIndex + (iStepInd * slideWindowStep) + epochOffset) * ones(nTrial, 1);
      iEpochEnd = iEpochBegin + slideWindowWidth;
      epochDuration = iEpochEnd - iEpochBegin;
      
      
      if iEpochEnd(1) > length(alignedRasters{1})
          break
      end
      switch dataType
         case 'spikes'
            nSpike = cellfun(@(x,y,z) sum(x(y:z)), alignedRasters, num2cell(iEpochBegin), num2cell(iEpochEnd), 'uniformoutput', false);
            iSpikeRate = cell2mat(nSpike) .* 1000 ./ epochDuration;
            % Choice dependence
            leftMetric = iSpikeRate(leftTrial);
            rightMetric = iSpikeRate(rightTrial);
         case 'eeg'
            eegMeanEpoch = cellfun(@(x,y,z) nanmean(x(y:z)), num2cell(alignedSignal,2), num2cell(iEpochBegin), num2cell(iEpochEnd), 'uniformoutput', false);
            eegMeanEpoch = cell2mat(eegMeanEpoch);
            % Choice dependence
            leftMetric = eegMeanEpoch(leftTrial);
            rightMetric = eegMeanEpoch(rightTrial);
         otherwise
      end
      
      % Break out of the loop if we're so far out that tChoice is
      % meaningless
      if (sum(isnan(leftMetric)) / length(leftMetric) > nanThreshold || ...
            sum(isnan(rightMetric)) / length(rightMetric) > nanThreshold)
         break
      end
      
      
      [p, h, stats] = ranksum(leftMetric, rightMetric);
      
      if p < alphaChoice
         choiceDependenceFound = true;
         tChoice = (iStepInd * slideWindowStep) + epochOffset;
      end
      iStepInd = iStepInd + 1;
   end
end



ddmData.choiceDependent     = choiceDependent;
ddmData.coherenceDependent  = coherenceDependent;
ddmData.ddmLike             = ddmLike;
ddmData.tChoice              = tChoice;
ddmData.leftIsIn            = leftIsIn;
ddmData.coeffIn                 = coeffIn;
ddmData.coeffOut                = coeffOut;
