function O = teg_repeated_measures_ANOVA(varargin)

% function O = teg_repeated_measures_ANOVA(M, levels, varnames, Betw_group, Betw_labels[, plots, fname])
%
% M is observation x variable-combination matrix.
% levels is vector of levels per factor.
% varnames is a cell array of strings.
%
% Other inputs are currently not intended for serious use!
%
% Maximum of 10 between-group levels.
%
% Thomas E. Gladwin (2010).

betw_subj_interaction_depth = 1;
pCritForFurther = 1;
% fprintf(['pCrit = ' num2str(pCritForFurther) '\n']);
only_show_FDR = 0;
verbose0 = 1;
plots0 = [];
digprec = 3; % digital precision of printouts.
doMANOVA = 0;

M = varargin{1};
levels = varargin{2};
p_fwe = 0.05 / prod(levels);
if size(M, 2) > prod(levels), % Number of observations per cell added
    b = size(M, 2) / 2;
    NM = M(:, (b + 1):end);
    M = M(:, 1:b);
else
    NM = ones(size(M));
end;
[nSubj, nVar] = size(M);
% Remove and store subject effects.
% Take different cell counts into account here.
subjEffect = [];
for iSubj = 1:size(M, 1),
    fNonNaN = ~isnan(M(iSubj, :));
    subjEffect(iSubj, 1) = sum(M(iSubj, fNonNaN) .* NM(iSubj, fNonNaN), 2) ./ sum(NM(iSubj, fNonNaN));
end;
M_subj = subjEffect * ones(1, size(M, 2));
M_raw = M;
M = M - M_subj;

varnames = varargin{3};
if isempty(varnames),
    varnames = cellstr(num2str((1:length(levels))'));
end;

groupVar = varargin{4};
groupVarnames = varargin{5};

if length(varargin) > 6,
    fname = varargin{7};
    if exist([fname '.ps']),
        delete([fname '.ps']);
    end;
else
    fname = [];
end;

if length(varargin) > 5,
    plots = varargin{6};
else
    plots = 0;
end;
if ~isempty(plots0),
    plots = plots0;
end;

% Get dummy-style matrices.
[X1, factorStarts, nColsFactor, labels, cellsets, factors] = teg_create_ANOVA_dummy(levels, varnames);
y = [];
y_raw = [];
yN = [];
X = [];
for iVar = 1:nVar,
    y = [y; M(:, iVar)];
    y_raw = [y_raw; M_raw(:, iVar)];
    yN = [yN; NM(:, iVar)];
    Xsub = ones(nSubj, 1) * X1(iVar, :);
    X = [X; Xsub];
end;

O.R = [];
O.labels = {};
curr_test = [];

for iPred = 1:length(factorStarts),
    predcols = factorStarts(iPred):(factorStarts(iPred) + nColsFactor(iPred) - 1);
    X0 = X(:, predcols);
    X_red_M = [];
    for iX = 1:size(X0, 2),
        [X_red_M0, dum] = teg_inner_recode_raw(reshape(X0(:, iX), size(M_raw)), ones(size(M_raw)), cellsets, iPred);
        X_red_M = [X_red_M X_red_M0(1, :)'];
    end;
    
    [y_red_M, dum] = teg_inner_recode_raw(M_raw, NM, cellsets, iPred);
    M_corr = M_raw;
    subjM = teg_mean(M_corr, 2);
    M_corr = M_corr - subjM * ones(1, size(M_corr, 2));
    [dum, expanded0] = teg_inner_recode_raw(M_corr, NM, cellsets, iPred);
    
    % Between-group
    nBetwInt = 2 ^ (size(groupVar, 2) - 1);
    for iBetwInt = 0:nBetwInt,
        coderStr = dec2bin(iBetwInt);
        betwStr = '';
        if iBetwInt == 0,
            BGv = [];
            DFB = 0;
        else,
            coder = (coderStr == '1');
            Bcols = find(coder);
            betwGroupThis = groupVar(:, Bcols);
            nLevels = [];
            combiner = [];
            for iBtmp = 1:size(betwGroupThis, 2),
                u = unique(betwGroupThis(:, iBtmp));
                nLevels = [nLevels length(u)];
                combiner = [combiner 10 .^ (iBtmp - 1)];
                betwStr = [betwStr ' x ' groupVarnames{Bcols(iBtmp)}];
            end;
            BGv = betwGroupThis * combiner(:);
            DFB = prod(nLevels - 1);
        end;
        [F, df1, df2, p, SSM, SSE, MSM, MSE, eta, eps0, bModel, X_ind] = teg_RMF_ANOVA(expanded0, X0, X0, BGv, DFB);
        O.R = [O.R; F df1 df2 p MSM MSE];
        thisLabel = [labels{iPred} betwStr];
        O.labels{length(O.labels) + 1} = thisLabel;
        
        if p <= pCritForFurther,
            teg_report('', labels, iPred, df1, digprec, df2, F, p, SSM, SSE, MSM, MSE, eta, eps0, p_fwe, y_red_M, fname, verbose0, bModel, X_red_M, pCritForFurther, thisLabel);
        end;
    end;
end;
