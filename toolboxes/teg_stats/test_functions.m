function O = test_functions

nSubj = 50;
levels = [2 2 3];
profile0 = zeros(1, prod(levels));
for n = 1:length(levels),
    varnames{n} = ['Var' num2str(n)];
end;
Betw_cont = randn(nSubj, 1);
Betw_labels{1} = 'B1';

profile0(1:levels(end):end) = 3;
M = randn(nSubj, prod(levels)) + ones(nSubj, 1) * profile0;
tmp = M(:, 1:levels(end):end);
M(:, 1:levels(end):end) = tmp .* (Betw_cont * ones(1, size(tmp, 2)));

M = M + randn(size(M));

O = teg_repeated_measures_ANOVA(M, levels, varnames, Betw_cont, Betw_labels);
