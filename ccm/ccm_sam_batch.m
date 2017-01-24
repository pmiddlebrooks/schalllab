%%

% modelType = 'noise1GoVary';
% modelType = 'noiseFreeGoVary';

SAM.des.noiseBound = true;
SAM.des.tBound = false;

if ~SAM.des.noiseBound && ~SAM.des.tBound
    modelType = 'noiseFreeTFree';
elseif SAM.des.noiseBound && ~SAM.des.tBound
    modelType = 'noiseBoundTFree';
elseif ~SAM.des.noiseBound && SAM.des.tBound
    modelType = 'noiseFreeTBound';
elseif SAM.des.noiseBound && SAM.des.tBound
    modelType = 'noiseBoundTBound';
end
%%

iSubj = 'broca';
SAM.des.condParam = 't0';
SAM.des.choiceMech.type = 'race';
goOrAll = 'go';
sam_spec_job

%%
ccm_sam_simulation

%%
localFigurePath = local_figure_path;
localFigurePath = [localFigurePath, 'sam/', modelType, '/'];
if ~isdir(localFigurePath)
    mkdir(localFigurePath)
end

print(55,[localFigurePath, iSubj, '_',goOrAll,'_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam,'_CDF_Def'],'-dpdf', '-r300')
print(57,[localFigurePath, iSubj, '_',goOrAll,'_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_CDF_Full'],'-dpdf', '-r300')
print(56,[localFigurePath, iSubj, '_',goOrAll,'_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_MeanRT'],'-dpdf', '-r300')


               %%         
        simple_plot_script_dynamics
        
                        figureHandle = 650;
        print(650,[localFigurePath, iSubj, '_',goOrAll, '_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_Traject_',mat2str(iCnd),'_',mat2str(iStopTrType)],'-dpdf', '-r300')



%%
localFigurePath = local_figure_path;
localFigurePath = [localFigurePath, 'sam/', modelType, '/'];
if ~isdir(localFigurePath)
    mkdir(localFigurePath)
end

print(55,[localFigurePath, iSubj, '_',goOrAll,'_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam,'_CDF_Def'],'-dpdf', '-r300')
print(57,[localFigurePath, iSubj, '_',goOrAll,'_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_CDF_Full'],'-dpdf', '-r300')
print(56,[localFigurePath, iSubj, '_',goOrAll,'_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_MeanRT'],'-dpdf', '-r300')


%%
                        figureHandle = 650;
simple_plot_script_dynamics
        print(650,[localFigurePath, iSubj, '_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_Traject_',mat2str(iCnd),'_',mat2str(iStopTrType)],'-dpdf', '-r300')

        %%
localFigurePath = local_figure_path;
localFigurePath = [localFigurePath, 'sam/', modelType, '/'];
if ~isdir(localFigurePath)
    mkdir(localFigurePath)
end


goOrAll = 'all';

print(55,[localFigurePath, iSubj, '_',goOrAll,'_', SAM.des.choiceMech.type, '_p',SAM.des.condParam,'_CDF'],'-dpdf', '-r300')
print(57,[localFigurePath, iSubj, '_',goOrAll,'_', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_CDF_Full'],'-dpdf', '-r300')
print(56,[localFigurePath, iSubj, '_',goOrAll,'_', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_MeanRT'],'-dpdf', '-r300')

%%
for iCnd = 1 : nSignal
    
    for iStopTrType = 2 : 4
        simple_plot_script_dynamics
        print(650,[localFigurePath, iSubj, '_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_Traject_',mat2str(iCnd),'_',mat2str(iStopTrType)],'-dpdf', '-r300')
        
        nTrialSim = 5;
        simple_plot_script_nTrial_dynamics
        print(651,[localFigurePath, iSubj, '_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_Traject_',mat2str(iCnd),'_',mat2str(iStopTrType),'_5Trial'],'-dpdf', '-r300')
        
        nTrialSim = 20;
        simple_plot_script_nTrial_dynamics
        print(651,[localFigurePath, iSubj, '_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_Traject_',mat2str(iCnd),'_',mat2str(iStopTrType),'_20Trial'],'-dpdf', '-r300')
        
        
        
    end
end

%%  START HERE
goOrAll = 'all';
paramMat = {'t0','v','zc'};
choiceMat = {'race','ffi','li'};
subMat = {'broca'};
modelConst ={'noiseFreeTFree','noiseBoundTFree','noiseFreeTBound','noiseBoundTBound'};
modelConst ={'noiseBoundTBound'};


for s = 1 : length(subMat)
    iSubj = subMat{s}
    
    for m = 1 : length(modelConst)
        modelConst{m}
        modelType = modelConst{m}
        switch modelConst{m};
            case 'noiseFreeTFree'
                SAM.des.noiseBound = false;
                SAM.des.tBound = false;
            case 'noiseBoundTFree'
                SAM.des.noiseBound = true;
                SAM.des.tBound = false;
            case 'noiseFreeTBound'
                SAM.des.noiseBound = false;
                SAM.des.tBound = true;
            case 'noiseBoundTBound'
                SAM.des.noiseBound = true;
                SAM.des.tBound = true;
        end
        
        for p = 1 : length(paramMat)
            paramMat{p}
            SAM.des.condParam = paramMat{p};
            
            for c = 1 : length(choiceMat)
                choiceMat{c}
                SAM.des.choiceMech.type = choiceMat{c};
                
                
                
                sam_spec_job
                
                ccm_sam_simulation
                
                
                localFigurePath = local_figure_path;
                localFigurePath = [localFigurePath, 'sam/', modelType, '/'];
                if ~isdir(localFigurePath)
                    mkdir(localFigurePath)
                end
                
                print(55,[localFigurePath, iSubj, '_',goOrAll,'_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam,'_CDF_Def'],'-dpdf', '-r300')
                print(57,[localFigurePath, iSubj, '_',goOrAll,'_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_CDF_Full'],'-dpdf', '-r300')
                print(56,[localFigurePath, iSubj, '_',goOrAll,'_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_MeanRT'],'-dpdf', '-r300')
                
                
                for iCnd = 3
                    
                    for iStopTrType = 3
                        figureHandle = 650;
                        simple_plot_script_dynamics
                        print(figureHandle,[localFigurePath, iSubj, '_',goOrAll, '_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_Traject_',mat2str(iCnd),'_',mat2str(iStopTrType)],'-dpdf', '-r300')
                        
                        figureHandle = 651;
                        nTrialSim = 5;
                        simple_plot_script_nTrial_dynamics
                        print(figureHandle,[localFigurePath, iSubj, '_',goOrAll, '_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_Traject_',mat2str(iCnd),'_',mat2str(iStopTrType),'_5Trial'],'-dpdf', '-r300')
                        
                        figureHandle = 652;
                        nTrialSim = 20;
                        simple_plot_script_nTrial_dynamics
                        print(figureHandle,[localFigurePath, iSubj, '_',goOrAll, '_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_Traject_',mat2str(iCnd),'_',mat2str(iStopTrType),'_20Trial'],'-dpdf', '-r300')
                        
                        
                        
                    end
                end
                
                
            end
        end
    end
end

goOrAll = 'all';
paramMat = {'t0','v','zc'};
choiceMat = {'race','ffi','li'};
subMat = {'xena', 'human'};
modelConst ={'noiseFreeTFree','noiseBoundTFree','noiseFreeTBound','noiseBoundTBound'};
modelConst ={'noiseBoundTFree','noiseFreeTBound','noiseBoundTBound'};

for s = 1 : length(subMat)
    iSubj = subMat{s}
    
    for m = 1 : length(modelConst)
        modelConst{m}
        modelType = modelConst{m}
        switch modelConst{m};
            case 'noiseFreeTFree'
                SAM.des.noiseBound = false;
                SAM.des.tBound = false;
            case 'noiseBoundTFree'
                SAM.des.noiseBound = true;
                SAM.des.tBound = false;
            case 'noiseFreeTBound'
                SAM.des.noiseBound = false;
                SAM.des.tBound = true;
            case 'noiseBoundTBound'
                SAM.des.noiseBound = true;
                SAM.des.tBound = true;
        end
        
        for p = 1 : length(paramMat)
            paramMat{p}
            SAM.des.condParam = paramMat{p};
            
            for c = 1 : length(choiceMat)
                choiceMat{c}
                SAM.des.choiceMech.type = choiceMat{c};
                
                
                
                sam_spec_job
                
                ccm_sam_simulation
                
                
                localFigurePath = local_figure_path;
                localFigurePath = [localFigurePath, 'sam/', modelType, '/'];
                if ~isdir(localFigurePath)
                    mkdir(localFigurePath)
                end
                
                print(55,[localFigurePath, iSubj, '_',goOrAll,'_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam,'_CDF_Def'],'-dpdf', '-r300')
                print(57,[localFigurePath, iSubj, '_',goOrAll,'_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_CDF_Full'],'-dpdf', '-r300')
                print(56,[localFigurePath, iSubj, '_',goOrAll,'_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_MeanRT'],'-dpdf', '-r300')
                
                
                for iCnd = 3
                    
                    for iStopTrType = 3
                        figureHandle = 650;
                        simple_plot_script_dynamics
                        print(figureHandle,[localFigurePath, iSubj, '_',goOrAll, '_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_Traject_',mat2str(iCnd),'_',mat2str(iStopTrType)],'-dpdf', '-r300')
                        
                        figureHandle = 651;
                        nTrialSim = 5;
                        simple_plot_script_nTrial_dynamics
                        print(figureHandle,[localFigurePath, iSubj, '_',goOrAll, '_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_Traject_',mat2str(iCnd),'_',mat2str(iStopTrType),'_5Trial'],'-dpdf', '-r300')
                        
                        figureHandle = 652;
                        nTrialSim = 20;
                        simple_plot_script_nTrial_dynamics
                        print(figureHandle,[localFigurePath, iSubj, '_',goOrAll, '_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_Traject_',mat2str(iCnd),'_',mat2str(iStopTrType),'_20Trial'],'-dpdf', '-r300')
                        
                        
                        
                    end
                end
                
                
            end
        end
    end
end



%%

goOrAll = 'all';
paramMat = {'t0','v','zc'};
choiceMat = {'race','ffi','li'};
subMat = {'xena', 'human'};
modelConst ={'noiseFreeTFree','noiseBoundTFree','noiseFreeTBound','noiseBoundTBound'};
modelConst ={'noiseFreeTFree'};

for s = 1 : length(subMat)
    
    for m = 1 : length(modelConst)
        modelType = modelConst{m}
        switch modelConst{m};
            case 'noiseFreeTFree'
                SAM.des.noiseBound = false;
                SAM.des.tBound = false;
            case 'noiseBoundTFree'
                SAM.des.noiseBound = true;
                SAM.des.tBound = false;
            case 'noiseFreeTBound'
                SAM.des.noiseBound = false;
                SAM.des.tBound = true;
            case 'noiseBoundTBound'
                SAM.des.noiseBound = true;
                SAM.des.tBound = true;
        end
        
        for p = 1 : length(paramMat)
            SAM.des.condParam = paramMat{p};
            
            for c = 1 : length(choiceMat)
                iSubj = subMat{s}
                modelConst{m}
                paramMat{p}
                choiceMat{c}
                SAM.des.choiceMech.type = choiceMat{c};
                
                
                
                sam_spec_job
                
                ccm_sam_simulation
                
                
                localFigurePath = local_figure_path;
                localFigurePath = [localFigurePath, 'sam/', modelType, '/'];
                if ~isdir(localFigurePath)
                    mkdir(localFigurePath)
                end
                
                print(55,[localFigurePath, iSubj, '_',goOrAll,'_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam,'_CDF_Def'],'-dpdf', '-r300')
                print(57,[localFigurePath, iSubj, '_',goOrAll,'_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_CDF_Full'],'-dpdf', '-r300')
                print(56,[localFigurePath, iSubj, '_',goOrAll,'_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_MeanRT'],'-dpdf', '-r300')
                
                
                for iCnd = 3
                    
                    for iStopTrType = 3
                        figureHandle = 650;
                        simple_plot_script_dynamics
                        print(figureHandle,[localFigurePath, iSubj, '_',goOrAll, '_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_Traject_',mat2str(iCnd),'_',mat2str(iStopTrType)],'-dpdf', '-r300')
                        
                        figureHandle = 651;
                        nTrialSim = 5;
                        simple_plot_script_nTrial_dynamics
                        print(figureHandle,[localFigurePath, iSubj, '_',goOrAll, '_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_Traject_',mat2str(iCnd),'_',mat2str(iStopTrType),'_5Trial'],'-dpdf', '-r300')
                        
                        figureHandle = 652;
                        nTrialSim = 20;
                        simple_plot_script_nTrial_dynamics
                        print(figureHandle,[localFigurePath, iSubj, '_',goOrAll, '_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_Traject_',mat2str(iCnd),'_',mat2str(iStopTrType),'_20Trial'],'-dpdf', '-r300')
                        
                        
                        
                    end
                end
                
                
            end
        end
    end
end

%%
goOrAll = 'all';
paramMat = {'t0','v''zc'};
choiceMat = {'race','ffi','li'};
subMat = {'xena','human'};
modelConst ={'noiseFreeTFree','noiseBoundTFree','noiseFreeTBound','noiseBoundTBound'};

for s = 1 : length(subMat)
    
    for m = 1 : length(modelConst)
        modelType = modelConst{m}
        switch modelConst{m};
            case 'noiseFreeTFree'
                SAM.des.noiseBound = false;
                SAM.des.tBound = false;
            case 'noiseBoundTFree'
                SAM.des.noiseBound = true;
                SAM.des.tBound = false;
            case 'noiseFreeTBound'
                SAM.des.noiseBound = false;
                SAM.des.tBound = true;
            case 'noiseBoundTBound'
                SAM.des.noiseBound = true;
                SAM.des.tBound = true;
        end
        
        for p = 1 : length(paramMat)
            SAM.des.condParam = paramMat{p};
            
            for c = 1 : length(choiceMat)
                SAM.des.choiceMech.type = choiceMat{c};
                
                
                
                
                iSubj = subMat{s}
                modelConst{m}
                paramMat{p}
                choiceMat{c}
                
                
                
                
                sam_spec_job
                
                ccm_sam_simulation
                
                
                localFigurePath = local_figure_path;
                localFigurePath = [localFigurePath, 'sam/', modelType, '/'];
                if ~isdir(localFigurePath)
                    mkdir(localFigurePath)
                end
                
                print(55,[localFigurePath, iSubj, '_',goOrAll,'_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam,'_CDF_Def'],'-dpdf', '-r300')
                print(57,[localFigurePath, iSubj, '_',goOrAll,'_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_CDF_Full'],'-dpdf', '-r300')
                print(56,[localFigurePath, iSubj, '_',goOrAll,'_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_MeanRT'],'-dpdf', '-r300')
                
                
                for iCnd = 3
                    
                    for iStopTrType = 3
                        figureHandle = 650;
                        simple_plot_script_dynamics
                        print(figureHandle,[localFigurePath, iSubj, '_',goOrAll, '_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_Traject_',mat2str(iCnd),'_',mat2str(iStopTrType)],'-dpdf', '-r300')
                        
                        figureHandle = 651;
                        nTrialSim = 5;
                        simple_plot_script_nTrial_dynamics
                        print(figureHandle,[localFigurePath, iSubj, '_',goOrAll, '_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_Traject_',mat2str(iCnd),'_',mat2str(iStopTrType),'_5Trial'],'-dpdf', '-r300')
                        
                        figureHandle = 652;
                        nTrialSim = 20;
                        simple_plot_script_nTrial_dynamics
                        print(figureHandle,[localFigurePath, iSubj, '_',goOrAll, '_c', SAM.des.choiceMech.type, '_p',SAM.des.condParam, '_Traject_',mat2str(iCnd),'_',mat2str(iStopTrType),'_20Trial'],'-dpdf', '-r300')
                    end
                end
                
           


            end
        end
    end
end


%%
goOrAll = 'go';
paramMat = {'v','zc'};
choiceMat = {'race','ffi','li'};
subMat = {'broca'};
modelConst ={'noiseBoundTBound'};

for s = 1 : length(subMat)
    
    for m = 1 : length(modelConst)
        modelType = modelConst{m}
        switch modelConst{m};
            case 'noiseFreeTFree'
                SAM.des.noiseBound = false;
                SAM.des.tBound = false;
            case 'noiseBoundTFree'
                SAM.des.noiseBound = true;
                SAM.des.tBound = false;
            case 'noiseFreeTBound'
                SAM.des.noiseBound = false;
                SAM.des.tBound = true;
            case 'noiseBoundTBound'
                SAM.des.noiseBound = true;
                SAM.des.tBound = true;
        end
        
        for p = 1 : length(paramMat)
            SAM.des.condParam = paramMat{p};
            
            for c = 1 : length(choiceMat)
                SAM.des.choiceMech.type = choiceMat{c};
                
                
                
                
                iSubj = subMat{s}
                modelConst{m}
                paramMat{p}
                choiceMat{c}
                
                
                
                
                sam_spec_job
                


            end
        end
    end
end

goOrAll = 'go';
paramMat = {'t0','v','zc'};
choiceMat = {'race','ffi','li'};
subMat = {'xena','human'};
modelConst ={'noiseBoundTBound'};

for s = 1 : length(subMat)
    
    for m = 1 : length(modelConst)
        modelType = modelConst{m}
        switch modelConst{m};
            case 'noiseFreeTFree'
                SAM.des.noiseBound = false;
                SAM.des.tBound = false;
            case 'noiseBoundTFree'
                SAM.des.noiseBound = true;
                SAM.des.tBound = false;
            case 'noiseFreeTBound'
                SAM.des.noiseBound = false;
                SAM.des.tBound = true;
            case 'noiseBoundTBound'
                SAM.des.noiseBound = true;
                SAM.des.tBound = true;
        end
        
        for p = 1 : length(paramMat)
            SAM.des.condParam = paramMat{p};
            
            for c = 1 : length(choiceMat)
                SAM.des.choiceMech.type = choiceMat{c};
                
                
                
                
                iSubj = subMat{s}
                modelConst{m}
                paramMat{p}
                choiceMat{c}
                
                
                
                
                sam_spec_job
                


            end
        end
    end
end
