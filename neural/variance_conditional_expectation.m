function varCE = variance_conditional_expectation(raster)

nSpike = nansum(raster, 2);

totalVar = var(nSpike);


% Need to find minimum fano factor, using 60 ms time windows of spike data.
% But how did Churchland et al 2006 do this? Did they run through the whole
% trial ms by ms, repeating it for all possible alignment times?

% fanoFactor = fano_factor(raster)


phi = min(fanoFactor);
ppv = point_process_variance(nSpike, phi);

varCE = totalVar - ppv;


function ppv = point_process_variance(nSpike, phi)

ppv = phi .* mean(nSpike);
end

end