function [m, se] = teg_make_table(M, digprec1, digprec2)

[N, nCol] = size(M);
m = teg_nanmean(M);
se = sqrt(teg_nanvar(M)) ./ sqrt(N);

for iCol = 1:nCol,
    fprintf([num2str(m(iCol), digprec1) ' (' num2str(se(iCol), digprec2) ')']);
    if iCol < nCol,
        fprintf('\t');
    else
        fprintf('\n');
    end;
end;
