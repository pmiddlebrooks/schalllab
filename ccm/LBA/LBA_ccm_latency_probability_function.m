function predictedRT = LBA_ccm_latency_probability_function(correct, incorrect, t, signalStrengthArray, FIT_STOPS)

% nCondition = size(fieldnames(correct));
% conditionArray = fieldnames(correct)
nCondition = length(signalStrengthArray);
if FIT_STOPS
    nGoStop = 2;
else
    nGoStop = 1;
end

for j = 1 : nGoStop  % Go, Stop
    for i = 1 : nCondition
        midW = (correct{j, i}(end) - correct{j, i}(1)) / 2;
        midWIndex = find(correct{j, i} > midW, 1);
        predictedRT.correct(j, i) = t{j, i}(midWIndex);
        
        midL = (incorrect{j, i}(end) - incorrect{j, i}(1)) / 2;
        midLIndex = find(incorrect{j, i} > midL, 1);
        predictedRT.incorrect(j, i) = t{j, i}(midLIndex);
        
    end
end