function sam_plot(SAM,prd);

obs = SAM.optim.obs;

trialCat    = SAM.optim.obs.trialCat;
nTrialCat   = size(trialCat,1);

minBinSize  = SAM.optim.cost.stat.minBinSize;

nHorPanel = 4;

% Set up panel
% figure('units','normalized','outerposition',[0 0 1 1]);
p = panel();
p.pack(ceil(nTrialCat/nHorPanel),nHorPanel);

plotType = 'normal';


for iTrialCat = 1:nTrialCat
    
    if ~isempty(regexp(trialCat{iTrialCat},'goTrial.*', 'once'))
        
        % Numbers
        nTotalPrd                       = prd.nTotal(iTrialCat);
        nGoCCorrPrd                     = prd.nGoCCorr(iTrialCat);
        nGoCErrorPrd                    = prd.nGoCError(iTrialCat);
        
        nGoCCorrObs                     = obs.nGoCCorr(iTrialCat);
        nGoCErrorObs                    = obs.nGoCError(iTrialCat);
        
        % Trial probabilities
        pGoCCorrPrd                     = prd.pGoCCorr(iTrialCat);
        pGoCErrorPrd                    = prd.pGoCError(iTrialCat);
        
        pGoCCorrObs                     = obs.pGoCCorr(iTrialCat);
        pGoCErrorObs                    = obs.pGoCError(iTrialCat);
        
        % Reaction times
        switch lower(plotType)
            case 'normal'
                rtGoCCorrPrd                    = prd.rtGoCCorr{iTrialCat};
                rtGoCErrorPrd                   = prd.rtGoCError{iTrialCat};
            case 'defective'
                rtGoCCorrPrd                    = nan(nTotalPrd,1);
                rtGoCCorrPrd(1:nGoCCorrPrd)     = prd.rtGoCCorr{iTrialCat};
                rtGoCErrorPrd                   = nan(nTotalPrd,1);
                rtGoCErrorPrd(1:nGoCErrorPrd)   = prd.rtGoCError{iTrialCat};
        end
        
        rtQGoCCorrObs                   = obs.rtQGoCCorr{iTrialCat};
        rtQGoCErrorObs                  = obs.rtQGoCError{iTrialCat};
        
        % Probabilities
        switch lower(plotType)
            case 'normal'
                cumProbGoCCor          = obs.cumProbGoCCorr{iTrialCat};
                cumProbGoCError        = obs.cumProbGoCError{iTrialCat};
            case 'defective'
                cumProbGoCCor          = obs.cumProbDefectiveGoCCorr{iTrialCat};
                cumProbGoCError        = obs.cumProbDefectiveGoCError{iTrialCat};
        end
        
        % Plot
        [iColumn,iRow] = ind2sub([nHorPanel,ceil(nTrialCat/nHorPanel)],iTrialCat);
        p(iRow,iColumn).select();
        p(iRow,iColumn).hold('on');
        
        % Predictions as lines
        if ~isempty(rtQGoCCorrObs)
            plot(rtGoCCorrPrd,cmtb_edf(rtGoCCorrPrd(:),rtGoCCorrPrd(:)),'Color','k','LineStyle','-');
        end
        
        if ~isempty(rtQGoCErrorObs)
            plot(rtGoCErrorPrd,cmtb_edf(rtGoCErrorPrd(:),rtGoCErrorPrd(:)),'Color','k','LineStyle','--');
        end
        
        % Observations as circles
        if ~isempty(rtQGoCCorrObs)
            plot(rtQGoCCorrObs,cumProbGoCCor,'ko');
        end
        
        if ~isempty(rtQGoCErrorObs)
            plot(rtQGoCErrorObs,cumProbGoCError,'kd');
        end
        
        % Print trial probabilities
        fprintf(1,'GoCCorrObs = %.2f, GoCCorrPrd = %.2f \n',pGoCCorrObs,pGoCCorrPrd);
        fprintf(1,'GoCErrorObs = %.2f, GoCErrorPrd = %.2f \n',pGoCErrorObs,pGoCErrorPrd);
        
    elseif ~isempty(regexp(trialCat{iTrialCat},'stopTrial.*', 'once'))
        
        % Numbers
        nTotalPrd                       = prd.nTotal(iTrialCat);
        nStopICorrPrd                   = prd.nStopICorr(iTrialCat);
        nStopIErrorCCorrPrd             = prd.nStopIErrorCCorr(iTrialCat);
        nStopIErrorCErrorPrd            = prd.nStopIErrorCError(iTrialCat);
        
        nStopICorrObs                   = obs.nStopICorr(iTrialCat);
        nStopIErrorCCorrObs             = obs.nStopIErrorCCorr(iTrialCat);
        nStopIErrorCErrorObs            = obs.nStopIErrorCError(iTrialCat);
        
        % Trial probabilities
        pStopICorrPrd                   = prd.pStopICorr(iTrialCat);
        pStopIErrorCCorrPrd             = prd.pStopIErrorCCorr(iTrialCat);
        pStopIErrorCErrorPrd            = prd.pStopIErrorCError(iTrialCat);
        
        pStopICorrObs                   = obs.pStopICorr(iTrialCat);
        pStopIErrorCCorrObs             = obs.pStopIErrorCCorr(iTrialCat);
        pStopIErrorCErrorObs            = obs.pStopIErrorCError(iTrialCat);
        
        
        % Reaction times
        switch lower(plotType)
            case 'normal'
                rtStopICorrPrd                              = prd.rtStopICorr{iTrialCat};
                rtStopIErrorCCorrPrd                        = prd.rtStopIErrorCCorr{iTrialCat};
                rtStopIErrorCErrorPrd                       = prd.rtStopIErrorCError{iTrialCat};
            case 'defective'
                rtStopICorrPrd                                 = nan(nTotalPrd,1);
                rtStopICorrPrd(1:nStopICorrPrd)                = prd.rtStopICorr{iTrialCat};
                rtStopIErrorCCorrPrd                           = nan(nTotalPrd,1);
                rtStopIErrorCCorrPrd(1:nStopIErrorCCorrPrd)    = prd.rtStopIErrorCCorr{iTrialCat};
                rtStopIErrorCErrorPrd                          = nan(nTotalPrd,1);
                rtStopIErrorCErrorPrd(1:rtStopIErrorCErrorPrd) = prd.rtStopIErrorCError{iTrialCat};
        end
        
        rtQStopIErrorCCorrObs                = obs.rtQStopIErrorCCorr{iTrialCat};
        rtQStopIErrorCErrorObs               = obs.rtQStopIErrorCError{iTrialCat};
        
        % Probabilities
        switch lower(plotType)
            case 'normal'
                cumProbStopIErrorCCorr          = obs.cumProbStopIErrorCCorr{iTrialCat};
                cumProbStopIErrorCError        = obs.cumProbStopIErrorCError{iTrialCat};
            case 'defective'
                cumProbStopIErrorCCorr          = obs.cumProbDefectiveStopIErrorCCorr{iTrialCat};
                cumProbStopIErrorCError        = obs.cumProbDefectiveStopIErrorCError{iTrialCat};
        end
        
        % Plot
        
        [iColumn,iRow] = ind2sub([nHorPanel,ceil(nTrialCat/nHorPanel)],iTrialCat);
        p(iRow,iColumn).select();
        p(iRow,iColumn).hold('on');
        
        % Predictions as lines
        plot(rtStopICorrPrd,cmtb_edf(rtStopICorrPrd(:),rtStopICorrPrd(:)),'Color','k','LineStyle','-.');
        
        if ~isempty(rtQStopIErrorCCorrObs)
            plot(rtStopIErrorCCorrPrd,cmtb_edf(rtStopIErrorCCorrPrd(:),rtStopIErrorCCorrPrd(:)),'Color','k','LineStyle','-');
        end
        
        if ~isempty(rtQStopIErrorCErrorObs)
            plot(rtStopIErrorCErrorPrd,cmtb_edf(rtStopIErrorCErrorPrd(:),rtStopIErrorCErrorPrd(:)),'Color','k','LineStyle','--');
        end
        
        % Observations as circles
        if ~isempty(rtQStopIErrorCCorrObs)
            plot(rtQStopIErrorCCorrObs,cumProbStopIErrorCCorr,'ko');
        end
        
        if ~isempty(rtQStopIErrorCErrorObs)
            plot(rtQStopIErrorCErrorObs,cumProbStopIErrorCError,'kd');
        end
        
        
        % Print trial probabilities
        fprintf(1,'StopICorrObs = %.2f, StopICorrPrd = %.2f \n',pStopICorrObs,pStopICorrPrd);
        fprintf(1,'StopIErrorCCorrObs = %.2f, StopIErrorCCorrPrd = %.2f \n',pStopIErrorCCorrObs,pStopIErrorCCorrPrd);
        fprintf(1,'StopIErrorCErrorObs = %.2f, StopIErrorCErrorPrd = %.2f \n',pStopIErrorCErrorObs,pStopIErrorCErrorPrd);
        
    end
            

  
%   rtPrdCorr               = nan(nTotalPrd,1);
%   rtPrdCorr(1:nPrdCorr)   = prd.rtCorr{iTrialCat};
%   rtPrdError              = nan(nTotalPrd,1);
%   rtPrdError(1:nPrdError) = prd.rtError{iTrialCat};
%   
%   nObsCorr                = obs.nCorr(iTrialCat);
%   nObsError               = obs.nError(iTrialCat);
%   
%   rtQObsCorr              = obs.rtQCorr{iTrialCat};
%   rtQObsError             = obs.rtQError{iTrialCat};
%   pDefectiveObsCorr       = obs.pDefectiveCorr{iTrialCat};
%   pDefectiveObsError      = obs.pDefectiveError{iTrialCat};
%   
%   % Plot
%   [iColumn,iRow] = ind2sub([nHorPanel,ceil(nTrialCat/nHorPanel)],iTrialCat);
%   p(iRow,iColumn).select();
%   p(iRow,iColumn).hold('on');
%   
%   % Predictions as lines
%   plot(rtPrdCorr,cmtb_edf(rtPrdCorr,rtPrdCorr),'Color','k','LineStyle','-');
%   plot(rtPrdError,cmtb_edf(rtPrdError,rtPrdError),'Color','k','LineStyle','--');
%   
%   % Observations as circles
%   plot(rtQObsCorr,pDefectiveObsCorr,'ko');
%   plot(rtQObsError,pDefectiveObsError,'kd');
%   
  % Set title
  title(obs.trialCat{iTrialCat});
  set(gca,'PlotBoxAspectRatio',[1.61,1,1],'XLim',[0 1200],'YLim',[0 1]);
    % Legend
%   legend(sprintf('Corr(N_O=%d,N_P=%d)',nObsCorr,nPrdCorr),sprintf('Error(N_O=%d,N_P=%d)',nObsError,nPrdError));
  
end