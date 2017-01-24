function ccm_n2pc_ddm_population(subjectID, varargin)


% Set defaults
task = 'ccm';
plotFlag    = 1;
printPlot   = 0;
figureHandle = 6200;
alignEvent = 'checkerOn';
electrodeArray  = eeg_electrode_map(subjectID);
for i = 1 : 2 : length(varargin)
   switch varargin{i}
      case 'plotFlag'
         plotFlag = varargin{i+1};
      case 'printPlot'
         printPlot = varargin{i+1};
      case 'figureHandle'
         figureHandle = varargin{i+1};
      case 'alignEvent'
         alignEvent = varargin{i+1};
      case 'electrodeArray'
         electrodeArray = varargin{i+1};
      otherwise
   end
end




switch subjectID
   case 'Human'
      pSignalArray = [.35 .42 .46 .5 .54 .58 .65];
   case 'Broca'
      pSignalArray = [.41 .45 .48 .5 .52 .55 .59];
   case 'Xena'
      pSignalArray = [.35 .42 .47 .5 .53 .58 .65];
end
signalLeftP = pSignalArray(pSignalArray < .5);
signalRightP = pSignalArray(pSignalArray > .5);


[sessionArray, subjectIDArray] = task_session_array(subjectID, task);


plotRange = [-200 : 300];
plotAlign = -plotRange(1);
nElectrode      = length(electrodeArray);
nSession = length(sessionArray);
channelArray = [];
for i = 1 : nElectrodem
   channelArray = [channelArray, find(strcmp(electrodeArray{i}, electrodeArray))];
end




% nRow = 2;
% nColumn = 2;
% axLGo = 1;
% axRGo = 2;
% axLStop = 3;
% axRStop = 4;
% figureHandle = 6124;
% [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = big_figure(nRow, nColumn, figureHandle, 'save');
% clf
% % [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_figure(nRow, nColumn, 'portrait', figureHandle);
% ax(axLGo) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 1) yAxesPosition(1, 1) axisWidth axisHeight]);
% ax(axRGo) = axes('units', 'centimeters', 'position', [xAxesPosition(1, 2) yAxesPosition(1, 2) axisWidth axisHeight]);
% ax(axLStop) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 1) yAxesPosition(2, 1) axisWidth axisHeight]);
% ax(axRStop) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 2) yAxesPosition(2, 2) axisWidth axisHeight]);
%
%
% %Colors
% rightColor = [1 0 0];
% leftColor = [0 0 1];
% stopContraColor = [1 0 0];
% stopIpsiColor = [0 0 1];
% diffColor = [0 0 0];
%
% hold(ax(axLGo), 'on')
% hold(ax(axRGo), 'on')
% hold(ax(axLStop), 'on')
% hold(ax(axRStop), 'on')




choiceErpL = cell(nSession, nElectrode);
choiceErpR = cell(nSession, nElectrode);
coherenceErp = cell(nSession, nElectrode, length(pSignalArray));

choiceErpMeanL = cell(nElectrode);
choiceErpMeanR = cell(nElectrode);
coherenceErpMean = cell(nElectrode);
% rtR             = nan(nSession, 1);
% rtL             = nan(nSession, 1);
emptyArray = [];
for iSession = 1 : nSession
   iSessionID = sessionArray{iSession};
   iSessionID
   iData = ccm_n2pc_ddm_like(subjectID, iSessionID, 'electrodeArray', electrodeArray, 'plotFlag', 0);
   
   if isempty(iData)
      emptyArray = [emptyArray, iSession];
      continue
   end
   
   
   % Loop through each electrode and analyze it
   for iChannelInd = 1 : nElectrode;
      iChannel = channelArray(iChannelInd);
      
      alignedSignal = iData(iChannel).alignedSignal;
   % Assign the variables common to each electrode we're analyzing
   leftTrial   = iData(iChannel).leftTrial;
   rightTrial  = iData(iChannel).rightTrial;
   signalP     = iData(iChannel).signalP;
   alignIndex  = iData(iChannel).alignIndex;
   epochOffset = iData(iChannel).epochOffset;
      
      
      
      % CHOICE DEPENDENCE: all left vs. all right trials
      iLeftSignal = nanmean(alignedSignal(leftTrial, :), 1);
      iRightSignal = nanmean(alignedSignal(rightTrial, :), 1);
      
      choiceErpL{iSession, iChannel} = nanmean(alignedSignal(leftTrial, alignIndex + plotRange), 1);
      choiceErpR{iSession, iChannel} = nanmean(alignedSignal(rightTrial, alignIndex + plotRange), 1);
      
      
      
      % COHERENCE DEPENDENCE
      for j = 1 : length(pSignalArray)
         jProp = pSignalArray(j);
         signalTrial = signalP == jProp;
         coherenceErp{iSession, iChannel, j} = ...
            nanmean(alignedSignal(signalTrial, alignIndex + plotRange), 1);
      end % Coherence loop
   end
   
   
   
   %    rtR(iSession) = nanmean(iData.rtR);
   %    rtL(iSession) = nanmean(iData.rtL);
   
end



% Loop through each electrode again and collect the means (across
% sessions) of the ERPs, and plot the data
for iChannelInd = 1 : nElectrode;
   iChannel = channelArray(iChannelInd);
   
   % CHOICE DEPENDENCE: all left vs. all right trials
   choiceErpMeanL{iChannel} = nanmean(cell2mat(choiceErpL(:, iChannel)), 1);
   choiceErpMeanR{iChannel} = nanmean(cell2mat(choiceErpR(:, iChannel)), 1);
   
   
   % COHERENCE DEPENDENCE
   for j = 1 : length(pSignalArray)
      coherenceErpMean{iChannel, j} = nanmean(cell2mat(coherenceErp(:, iChannel, j)), 1);
   end
   
   
   
   
   
   % ________________________________________________________________
   % PLOT THE DATA
   
   if plotFlag
      
      % SET UP PLOT
      lineW = 2;
      cMap = ccm_colormap(pSignalArray);
      leftColor = cMap(1,:) .* .8;
      rightColor = cMap(end,:) .* .8;
      nRow = 2;
      nColumn = 3;
      figureHandle = figureHandle + 1;
      if printPlot
         [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_landscape(nRow, nColumn, figureHandle);
      else
         [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = screen_figure(nRow, nColumn, figureHandle);
      end
      clf
      axChoice = 1;
      axCoh = 2;
      axCohL = 3;
      axCohR = 4;
      
      
      
      ax(axChoice) = axes('units', 'centimeters', 'position', [xAxesPosition(axChoice, 2) yAxesPosition(axChoice, 2) axisWidth axisHeight]);
      cla
      hold(ax(axChoice), 'on')
%       switch choiceDependent(iChannel)
%          case true
%             choiceStr = 'YES';
%          otherwise
%             choiceStr = 'NO';
%       end
%       tt = sprintf('Choice dependence: %s', choiceStr);
%       title(tt)
      
      ax(axCohL) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 1) yAxesPosition(2, 1) axisWidth axisHeight]);
      cla
      hold(ax(axCohL), 'on')
      title('Coherence dependence')
      
      ax(axCohR) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 3) yAxesPosition(2, 3) axisWidth axisHeight]);
      cla
      hold(ax(axCohR), 'on')
      title('Coherence dependence')
      
      ax(axCoh) = axes('units', 'centimeters', 'position', [xAxesPosition(2, 2) yAxesPosition(2, 2) axisWidth axisHeight]);
      cla
      hold(ax(axCoh), 'on')
%       switch coherenceDependent(iChannel)
%          case true
%             cohStr = 'YES';
%          otherwise
%             cohStr = 'NO';
%       end
%       tt = sprintf('Coherence dependence: %s', cohStr);
%       title(tt)
      
      
      signalMax = max(max(choiceErpMeanL{iChannel}), max(choiceErpMeanL{iChannel}));
      yMax = 1.1 * signalMax;
      signalMin = min(min(choiceErpMeanL{iChannel}), min(choiceErpMeanL{iChannel}));
      yMin = 1.1 * signalMin;
      fillX = [epochOffset, epochOffset +  100, epochOffset + 100, epochOffset];
      fillY = [yMin yMin yMax yMax];
      fillColor = [1 1 .5];
      
      
      
      
      
      % CHOICE DEPENDENCE PLOTTING(LEFT VS. RIGHT CHOICE FOR CORRECT TRIALS)
      axes(ax(axChoice))
      h = fill(fillX, fillY, fillColor);
      set(h, 'edgecolor', 'none');
      plot(ax(axChoice), plotRange, choiceErpMeanL{iChannel}, 'color', leftColor, 'linewidth', lineW)
      plot(ax(axChoice), plotRange, choiceErpMeanR{iChannel}, 'color', rightColor, 'linewidth', lineW)
      plot(ax(axChoice), [1 1], [yMin yMax], '-k', 'linewidth', 2);
      set(ax(axChoice), 'ylim', [yMin yMax])
      
      
      
      
      
      % COHERENCE DEPENDENCE PLOTTING
      
      % Leftward trials
      plot(ax(axCohL), [0 0], [yMin yMax], '-k', 'linewidth', 2);
      axes(ax(axCohL))
      h = fill(fillX, fillY, fillColor);
      set(h, 'edgecolor', 'none');
      
      for j = 1 : length(signalLeftP)
         
         % Determine color to use for plot based on which checkerboard color
         % proportion being used. Normalize the available color spectrum to do
         % it
         sigColor = cMap(j,:);
         
         plot(ax(axCohL), plotRange, coherenceErpMean{iChannel, j}, 'color', sigColor, 'linewidth', lineW)
         set(ax(axCohL), 'ylim', [yMin yMax])
         
%          scatter(ax(axCoh), trialDataSub.targ1CheckerProp(signalTrial), eegMeanEpoch(signalTrial), 'o', 'markeredgecolor', sigColor, 'markerfacecolor', sigColor, 'sizedata', 20)
      end % for i = 1 : length(signalLeftP)
      
      
      
      % Rightward trials
      plot(ax(axCohR), [0 0], [yMin yMax], '-k', 'linewidth', 2);
      axes(ax(axCohR))
      h = fill(fillX, fillY, fillColor);
      set(h, 'edgecolor', 'none');
      for j = (j+1) : (length(signalLeftP) + length(signalRightP))
         iProp = pSignalArray(j);
         
         % Determine color to use for plot based on which checkerboard color
         % proportion being used. Normalize the available color spectrum to do
         % it
         sigColor = cMap(j,:);
         
         plot(ax(axCohR), plotRange,coherenceErpMean{iChannel, j}, 'color', sigColor, 'linewidth', lineW)
         set(ax(axCohR), 'ylim', [yMin yMax])
         
%          scatter(ax(axCoh), trialDataSub.targ1CheckerProp(signalTrial), eegMeanEpoch(signalTrial), 'o', 'markeredgecolor', sigColor, 'markerfacecolor', sigColor, 'sizedata', 30)
      end % for i = 1 : length(signalRightP)
      
      
      
            h=axes('Position', [0 0 1 1], 'Visible', 'Off');
%       if choiceDependent(iChannel) && coherenceDependent(iChannel)
%          ddmStr = 'YES';
%       else
%          ddmStr = 'NO';
%       end
%       titleString = sprintf('%s \t %s \t DDM-Like: %s', sessionID, electrodeArray{iChannel}, ddmStr);
      titleString = sprintf('Across Sessions \t Electrode: %s ', electrodeArray{iChannel});
      text(0.5,1, titleString, 'HorizontalAlignment','Center', 'VerticalAlignment','Top')

      if printPlot
         localFigurePath = local_figure_path;
         print(figureHandle,[localFigurePath, subjectID, '_ccm_ddm_like_pop_', electrodeArray{iChannel}],'-dpdf', '-r300')
      end
   end % if plotFlag
end

% 
% 
% 
% 
% 
% 
% displayRange = iData.displayRange;
% alignTime = iData.alignTime;
% % Remove sessions that did not collect data;
% o1RightTargEEG(emptyArray) = [];
% o1LeftTargEEG(emptyArray) = [];
% o2LeftTargEEG(emptyArray) = [];
% o2RightTargEEG(emptyArray) = [];
% rtR(emptyArray) = [];
% rtL(emptyArray) = [];
% 
% lContraMean = nanmean(cell2mat(cellfun(@(x) nanmean(x, 1), o1RightTargEEG, 'uniformoutput', false)), 1);
% lIpsiMean = nanmean(cell2mat(cellfun(@(x) nanmean(x, 1), o1LeftTargEEG, 'uniformoutput', false)), 1);
% olDiff = lContraMean - lIpsiMean;
% 
% rContraMean = nanmean(cell2mat(cellfun(@(x) nanmean(x, 1), o2LeftTargEEG, 'uniformoutput', false)), 1);
% rIpsiMean = nanmean(cell2mat(cellfun(@(x) nanmean(x, 1), o2RightTargEEG, 'uniformoutput', false)), 1);
% o2Diff = rContraMean - rIpsiMean;
% 
% 
% rtMeanR = mean(rtR);
% rtMeanL = mean(rtL);
% % lContraMean = mean(cell2mat(o1RightTargEEG), 1);
% % lIpsiMean = mean(cell2mat(o1LeftTargEEG), 1);
% % olDiff = lContraMean - lIpsiMean;
% %
% % rContraMean = mean(cell2mat(o2LeftTargEEG), 1);
% % rIpsiMean = mean(cell2mat(o2RightTargEEG), 1);
% % orDiff = rContraMean - rIpsiMean;
% 
% % xMax = 1.1 * max(rtMeanL, rtMeanR);
% xMax = 500;
% 
% if plotFlag
%    plot(ax(axLGo), lContraMean, 'color', rightColor, 'linewidth', 2)
%    plot(ax(axLGo), lIpsiMean, 'color', leftColor, 'linewidth', 2)
%    plot(ax(axLGo), olDiff, 'color', diffColor, 'linewidth', 2)
%    ylim(ax(axLGo), [-.015 .015])
%    plot(ax(axLGo), [alignTime alignTime], ylim(ax(axLGo)))
%    xlim(ax(axLGo), [0 xMax])
%    plot(ax(axLGo), [rtMeanR + -displayRange(1)  rtMeanR + -displayRange(1)], ylim(ax(axLGo)), '--', 'color', rightColor)
%    plot(ax(axLGo), [rtMeanL + -displayRange(1)  rtMeanL + -displayRange(1)], ylim(ax(axLGo)), '--', 'color', leftColor)
%    set(ax(axLGo), 'XTickLabel', (-alignTime : 100 : xMax))
%    legend(ax(axLGo), 'o1 Right (Contra)', 'o1 Left (Ipsi)', 'Contra - Ipsi')
%    
%    plot(ax(axRGo), rContraMean, 'color', leftColor, 'linewidth', 2)
%    plot(ax(axRGo), rIpsiMean, 'color', rightColor, 'linewidth', 2)
%    plot(ax(axRGo), o2Diff, 'color', diffColor, 'linewidth', 2)
%    ylim(ax(axRGo), [-.015 .015])
%    plot(ax(axRGo), [alignTime alignTime], ylim(ax(axRGo)))
%    xlim(ax(axRGo), [0 xMax])
%    plot(ax(axRGo), [rtMeanR + -displayRange(1)  rtMeanR + -displayRange(1)], ylim(ax(axRGo)), '--', 'color', rightColor)
%    plot(ax(axRGo), [rtMeanL + -displayRange(1)  rtMeanL + -displayRange(1)], ylim(ax(axRGo)), '--', 'color', leftColor)
%    set(ax(axRGo), 'XTickLabel', (-alignTime : 100 : xMax))
%    legend(ax(axRGo), 'o2 Left (Contra)', 'o2 Right (Ipsi)', 'Contra - Ipsi')
% end
% 
% 
% data.o1RightTargEEG = o1RightTargEEG;
% data.o1LeftTargEEG = o1LeftTargEEG;
% data.o2LeftTargEEG = o2LeftTargEEG;
% data.o2RightTargEEG = o2RightTargEEG;
% 

