function sam_plot_obs_prd_inh_rt(SAM,prd);

defective = true;

obs = SAM.optim.obs;

nCnd  = SAM.expt.nCnd;
nSsd  = SAM.expt.nSsd;

nRows = nCnd;

if defective
    nCols = nSsd + 2; % Go trials +  Inhibition function + nSsd * signal-respond RTs
else
    nCols = 3;
end

iGoTrial = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,'^goTrial')),obs.trialCat,'Uni',0)));
iStopTrial = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,'^stopTrial')),obs.trialCat,'Uni',0)));

iStopTrial = reshape(iStopTrial,nSsd,nCnd)';

% Set up panel
% figure('units','normalized','outerposition',[0 0 1 1]);
p = panel('no-manage-font');
p.pack(nRows,nCols);

lnColor = {repmat([0 128 0],5,1), ...
           repmat([0 0 255],5,1), ...
           repmat([255 0 0],5,1)};

% lnColor = {[224 255 224; ...
%             168 224 168; ...
%             112 192 112; ...
%              56 160  56; ...
%               0 128   0], ...
%            [224 224 255; ...
%             168 168 255; ...
%             112 112 255; ...
%              56  56 255; ...
%               0   0 255], ...
%            [255 224 224; ...
%             255 192 192; ...
%             255 112 112; ...
%             255  56  56; ...
%             255   0   0]};
lnColor = cellfun(@(in1) in1./255,lnColor,'Uni',0);
lnStyleCorr = '-';
lnStyleError = '--';
lnWidth = 2;

for iCnd = 1:nCnd
    
    % Go RT data
    % ==========
    iDSRow = iGoTrial(iCnd);
    
    % Get the data
    [rtQObsCorr,rtQObsError,FObsQCorr,FObsQError, ...
     rtPrdCorr,rtPrdError,FPrdCorr,FPrdError] = get_data(defective,obs,prd,iDSRow);

    % Plot
    p(iCnd,1).select();
    p(iCnd,1).hold('on');
    
    % Predictions as lines
    hL1 = plot(rtPrdCorr,FPrdCorr,'Color',lnColor{iCnd}(end,:),'LineStyle',lnStyleCorr,'LineWidth',lnWidth);
    hL2 = plot(rtPrdError,FPrdError,'Color',lnColor{iCnd}(end,:),'LineStyle',lnStyleError,'LineWidth',lnWidth);

    % Observations as circles
    plot(rtQObsCorr,FObsQCorr,'ko','Color',lnColor{iCnd}(end,:));
    plot(rtQObsError,FObsQError,'kd','Color',lnColor{iCnd}(end,:));
    
    % Set axes
    set(gca,'PlotBoxAspectRatio',[1.61,1,1], ...
            'XLim',[0 2000], ...
            'XTick',[], ...
            'YLim',[0 1], ...
            'YTick',0:0.5:1, ...
            'YTickLabel','0.0|0.5|1.0');
    if iCnd == 1
        title('Go');
    end
    
    ylabel('Cum.prob.');
    
    p.fixdash(hL2,'--');
    
    if iCnd == 3
        xlabel('Time (ms)');
        set(gca,'PlotBoxAspectRatio',[1.61,1,1], ...
                'XTick',0:1000:2000);
    end
    
    % Inhibition function
    % ===================
    iDSRow = iStopTrial(iCnd,:);
    
    % Plot
    p(iCnd,2).select();
    p(iCnd,2).hold('on');
    
    % Predictions as lines
    plot(prd.ssd(iStopTrial(iCnd,:)),prd.pError(iStopTrial(iCnd,:)),'Color',lnColor{iCnd}(end,:),'LineStyle',lnStyleCorr,'LineWidth',lnWidth);
    
    % Observations as circles
    plot(obs.ssd(iStopTrial(iCnd,:)),obs.pError(iStopTrial(iCnd,:)),'ko','Color',lnColor{iCnd}(end,:));
    
    % Set axes
    set(gca,'PlotBoxAspectRatio',[1.61,1,1],'XLim',[0 1000],'YLim',[0 1]);
    if iCnd == 1
        title('Inhibition function');
    end
    
    for iSsd = 1:nSsd
        
        % Stop trial data
        % ===============
        iDSRow = iStopTrial(iCnd,iSsd);
        
         % Get the data
         [rtQObsCorr,rtQObsError,FObsQCorr,FObsQError, ...
          rtPrdCorr,rtPrdError,FPrdCorr,FPrdError] = get_data(defective,obs,prd,iDSRow);

        % Plot
        if defective
            p(iCnd,2 + iSsd).select();
            p(iCnd,2 + iSsd).hold('on');
        else
            p(iCnd,3).select();
            p(iCnd,3).hold('on');
        end
        
        % Predictions as lines
        hL1 = plot(rtPrdCorr,FPrdCorr,'Color',lnColor{iCnd}(iSsd,:),'LineStyle',lnStyleCorr,'LineWidth',lnWidth);
        hL2 = plot(rtPrdError,FPrdError,'Color',lnColor{iCnd}(iSsd,:),'LineStyle',lnStyleError,'LineWidth',lnWidth);

        % Observations as circles
        plot(rtQObsCorr,FObsQCorr,'ko','Color',lnColor{iCnd}(iSsd,:));
        plot(rtQObsError,FObsQError,'kd','Color',lnColor{iCnd}(iSsd,:));
        
        % Set axes
        set(gca,'PlotBoxAspectRatio',[1.61,1,1],'XLim',[0 2000],'YLim',[0 1]);
        if iCnd == 1
            title(sprintf('Stop, ssd %.0f ms',obs.ssd(iDSRow)));
            
        p.fixdash(hL2,'--');
        end
    end
end


function [rtQObsCorr,rtQObsError,FObsQCorr,FObsQError, ...
          rtPrdCorr,rtPrdError,FPrdCorr,FPrdError] = get_data(defective,obs,prd,iDSRow);

if defective
    
    nPrdTotal               = prd.nTotal(iDSRow);
    nPrdCorr                = prd.nCorr(iDSRow);
    nPrdError               = prd.nError(iDSRow);

    % Observed RTs
    rtQObsCorr              = obs.rtQCorr{iDSRow}(:);
    rtQObsError             = obs.rtQError{iDSRow}(:);
    
    % Observed cumulative probabilities
    FObsQCorr                = obs.pDefectiveCorr{iDSRow};
    FObsQError               = obs.pDefectiveError{iDSRow};
    
    % Predicted RTs
    rtPrdCorr               = nan(nPrdTotal,1);
    rtPrdCorr(1:nPrdCorr)   = prd.rtCorr{iDSRow}(:);
    rtPrdError              = nan(nPrdTotal,1);
    rtPrdError(1:nPrdError) = prd.rtError{iDSRow}(:);
    
    % Predicted cumulative probabilities
    FPrdCorr                = cmtb_edf(rtPrdCorr,rtPrdCorr);
    FPrdError               = cmtb_edf(rtPrdError,rtPrdError);
    
else
    
    % Observed RTs
    rtQObsCorr              = obs.rtQCorr{iDSRow}(:);
    rtQObsError             = obs.rtQError{iDSRow}(:);
    
    % Observed cumulative probabilities
    FObsQCorr                = obs.pDefectiveCorr{iDSRow}./obs.pCorr(iDSRow);
    FObsQError               = obs.pDefectiveError{iDSRow}./obs.pError(iDSRow);
    
    % Predicted RTs
    rtPrdCorr               = prd.rtCorr{iDSRow}(:);
    rtPrdError              = prd.rtError{iDSRow}(:);
    
    % Predicted cumulative probabilities
    FPrdCorr                = cmtb_edf(rtPrdCorr,rtPrdCorr);
    FPrdError               = cmtb_edf(rtPrdError,rtPrdError);
    
end