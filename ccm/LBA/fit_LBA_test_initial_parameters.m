function modelData = fit_LBA_test_initial_parameters(subjectID, sessionID)

AInitial = [100 200 300];
bInitial = [150 250 350];
vInitial = [.6 .7 .8];
T0Initial = [50 100];
sInitial = [.2 .3 .4];

ALL_FIXED_FLAG = 0;
FREE_BTWN_TARGET_FLAG = 2;
FREE_FLAG = 1;







% % Some Free
% A_freeOrFix = FREE_BTWN_TARGET_FLAG;
% b_freeOrFix = FREE_BTWN_TARGET_FLAG;
% v_freeOrFix = FREE_FLAG;
% t0_freeOrFix = FREE_BTWN_TARGET_FLAG;
% s_freeOrFix = ALL_FIXED_FLAG;

% All Free
A_freeOrFix = FREE_FLAG;
b_freeOrFix = FREE_FLAG;
v_freeOrFix = FREE_FLAG;
t0_freeOrFix = FREE_FLAG;
s_freeOrFix = ALL_FIXED_FLAG;


freeOrFixFlag = [A_freeOrFix, b_freeOrFix, v_freeOrFix, t0_freeOrFix, s_freeOrFix];

plotFlag = 1;







% initialParameter = {};
% solution = {};
% LL = {};
% AIC = {};
% BIC = {};
% CDF = {};
% correct = {};
% incorrect = {};
% t = {};
% 
% for AIndex = 1 : length(AInitial)
%     iA = AInitial(AIndex);
%     for bIndex = 1 : length(bInitial)
%         ib = bInitial(bIndex);
%         for vIndex = 1 : length(vInitial)
%             iv = vInitial(vIndex);
%             for T0Index = 1 : length(T0Initial)
%                 iT0 = T0Initial(T0Index);
%                 for sIndex = 1 : length(sInitial)
%                     is = sInitial(sIndex);
%                     
%                     fprintf('\nA: %.2f\t b: %.2f\t v: %.2f\t T0: %.2f\t s: %.2f\t\n', iA, ib, iv, iT0, is)
%                     InitialParamStruct.A = iA;
%                     InitialParamStruct.b = ib;
%                     InitialParamStruct.v = iv;
%                     InitialParamStruct.T0 = iT0;
%                     InitialParamStruct.s = is;
%                     [iSolution,iLL,iAIC,iBIC,iCDF,~, iCorrect, iIncorrect, iT] = fit_LBA_ccm(subjectID, sessionID, freeOrFixFlag, InitialParamStruct, plotFlag);
%                     initialParameter = [initialParameter; InitialParamStruct];
%                     solution = [solution; iSolution];
%                     LL = [LL; iLL];
%                     AIC = [AIC; iAIC];
%                     BIC = [BIC; iBIC];
%                     CDF = [CDF; iCDF];
%                     correct = [correct; iCorrect];
%                     incorrect = [incorrect; iIncorrect];
%                     t = [t; iT];
%                 end
%             end
%         end
%     end
% end
% size(initialParameter)
% size(solution)
% size(LL)
% size(AIC)
% size(BIC)
% size(CDF)
% size(correct)
% size(incorrect)
% size(t)
% 
% modelData = dataset(...
%     {initialParameter,                'initialParameter'},...
%     {solution,                'solution'},...
%     {LL,                'LL'},...
%     {AIC,                'AIC'},...
%     {BIC,                'BIC'},...
%     {CDF,                'CDF'},...
%     {correct,                'correct'},...
%     {incorrect,                'incorrect'},...
%     {t,  	't'});
% 
% save('test_initial_paramters_AllFreeLL.mat', 'modelData')





initialParameter = {};
solution = {};
chi2 = {};
cdfData = {};
cdfModel = {};
correct = {};
incorrect = {};
t = {};
for AIndex = 1 : length(AInitial)
    iA = AInitial(AIndex);
    for bIndex = find(bInitial >= iA)%1 : length(bInitial)
        ib = bInitial(bIndex);
        for vIndex = 1 : length(vInitial)
            iv = vInitial(vIndex);
            for T0Index = 1 : length(T0Initial)
                iT0 = T0Initial(T0Index);
                for sIndex = 1 : length(sInitial)
                    is = sInitial(sIndex);
                    
                    fprintf('\nA: %.2f\t b: %.2f\t v: %.2f\t T0: %.2f\t s: %.2f\t\n', iA, ib, iv, iT0, is)
                    InitialParamStruct.A = iA;
                    InitialParamStruct.b = ib;
                    InitialParamStruct.v = iv;
                    InitialParamStruct.T0 = iT0;
                    InitialParamStruct.s = is;
                    [iSolution,iChi2,iCdfData,iCdfModel, ~, iCorrect, iIncorrect, iT] = fit_LBA_chi2_ccm(subjectID, sessionID, freeOrFixFlag, InitialParamStruct, plotFlag);
                    initialParameter = [initialParameter; InitialParamStruct];
                    solution = [solution; iSolution];
                    chi2 = [chi2; iChi2];
                    cdfData = [cdfData; iCdfData];
                    cdfModel = [cdfModel; iCdfModel];
                    correct = [correct; iCorrect];
                    incorrect = [incorrect; iIncorrect];
                    t = [t; iT];
                end
            end
        end
    end
end

modelData = dataset(...
    {initialParameter,                'initialParameter'},...
    {solution,                'solution'},...
    {chi2,                'chi2'},...
    {cdfData,                'cdfData'},...
    {cdfModel,                'cdfModel'},...
    {correct,                'correct'},...
    {incorrect,                'incorrect'},...
    {t,  	't'});

save('test_initial_paramters_AllFreeChi_single.mat', 'modelData')




















% % SOME FREE
% A_freeOrFix = FREE_BTWN_TARGET_FLAG;
% b_freeOrFix = FREE_BTWN_TARGET_FLAG;
% v_freeOrFix = FREE_FLAG;
% t0_freeOrFix = FREE_BTWN_TARGET_FLAG;
% s_freeOrFix = ALL_FIXED_FLAG;
% 
% 
% freeOrFixFlag = [A_freeOrFix, b_freeOrFix, v_freeOrFix, t0_freeOrFix, s_freeOrFix];
% 
% plotFlag = 1;
% 
% 
% initialParameter = {};
% solution = {};
% chi2 = {};
% cdfData = {};
% cdfModel = {};
% correct = {};
% incorrect = {};
% t = {};
% for AIndex = 1 : length(AInitial)
%     iA = AInitial(AIndex);
%     for bIndex = 1 : length(bInitial)
%         ib = bInitial(bIndex);
%         for vIndex = 1 : length(vInitial)
%             iv = vInitial(vIndex);
%             for T0Index = 1 : length(T0Initial)
%                 iT0 = T0Initial(T0Index);
%                 for sIndex = 1 : length(sInitial)
%                     is = sInitial(sIndex);
%                     
%                     fprintf('\nA: %.2f\t b: %.2f\t v: %.2f\t T0: %.2f\t s: %.2f\t\n', iA, ib, iv, iT0, is)
%                     InitialParamStruct.A = iA;
%                     InitialParamStruct.b = ib;
%                     InitialParamStruct.v = iv;
%                     InitialParamStruct.T0 = iT0;
%                     InitialParamStruct.s = is;
%                     [iSolution,iChi2,iCdfData,iCdfModel, ~, iCorrect, iIncorrect, iT] = fit_LBA_chi2_ccm(subjectID, sessionID, freeOrFixFlag, InitialParamStruct, plotFlag);
%                     initialParameter = [initialParameter; InitialParamStruct];
%                     solution = [solution; iSolution];
%                     chi2 = [chi2; iChi2];
%                     cdfData = [cdfData; iCdfData];
%                     cdfModel = [cdfModel; iCdfModel];
%                     correct = [correct; iCorrect];
%                     incorrect = [incorrect; iIncorrect];
%                     t = [t; iT];
%                 end
%             end
%         end
%     end
% end
% 
% modelData = dataset(...
%     {initialParameter,                'initialParameter'},...
%     {solution,                'solution'},...
%     {chi2,                'chi2'},...
%     {cdfData,                'cdfData'},...
%     {cdfModel,                'cdfModel'},...
%     {correct,                'correct'},...
%     {incorrect,                'incorrect'},...
%     {t,  	't'});
% 
% save('test_initial_paramters_SomeFreeChi.mat', 'modelData')
% 
% 
% 
% initialParameter = {};
% solution = {};
% LL = {};
% AIC = {};
% BIC = {};
% CDF = {};
% correct = {};
% incorrect = {};
% t = {};
% 
% for AIndex = 1 : length(AInitial)
%     iA = AInitial(AIndex);
%     for bIndex = 1 : length(bInitial)
%         ib = bInitial(bIndex);
%         for vIndex = 1 : length(vInitial)
%             iv = vInitial(vIndex);
%             for T0Index = 1 : length(T0Initial)
%                 iT0 = T0Initial(T0Index);
%                 for sIndex = 1 : length(sInitial)
%                     is = sInitial(sIndex);
%                     
%                     fprintf('\nA: %.2f\t b: %.2f\t v: %.2f\t T0: %.2f\t s: %.2f\t\n', iA, ib, iv, iT0, is)
%                     InitialParamStruct.A = iA;
%                     InitialParamStruct.b = ib;
%                     InitialParamStruct.v = iv;
%                     InitialParamStruct.T0 = iT0;
%                     InitialParamStruct.s = is;
%                     [iSolution,iLL,iAIC,iBIC,iCDF,~, iCorrect, iIncorrect, iT] = fit_LBA_ccm(subjectID, sessionID, freeOrFixFlag, InitialParamStruct, plotFlag);
%                     initialParameter = [initialParameter; InitialParamStruct];
%                     solution = [solution; iSolution];
%                     LL = [LL; iLL];
%                     AIC = [AIC; iAIC];
%                     BIC = [BIC; iBIC];
%                     CDF = [CDF; iCDF];
%                     correct = [correct; iCorrect];
%                     incorrect = [incorrect; iIncorrect];
%                     t = [t; iT];
%                 end
%             end
%         end
%     end
% end
% size(initialParameter)
% size(solution)
% size(LL)
% size(AIC)
% size(BIC)
% size(CDF)
% size(correct)
% size(incorrect)
% size(t)
% 
% modelData = dataset(...
%     {initialParameter,                'initialParameter'},...
%     {solution,                'solution'},...
%     {LL,                'LL'},...
%     {AIC,                'AIC'},...
%     {BIC,                'BIC'},...
%     {CDF,                'CDF'},...
%     {correct,                'correct'},...
%     {incorrect,                'incorrect'},...
%     {t,  	't'});
% 
% save('test_initial_paramters_SomeFreeLL.mat', 'modelData')
