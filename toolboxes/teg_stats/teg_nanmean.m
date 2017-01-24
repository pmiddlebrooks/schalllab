function m = tegnanmean(M)

m = zeros(1, size(M, 2));
for iCol = 1:size(M, 2),
    vec = M(:, iCol);
    f = find(isnan(vec));
    vec(f) = [];
    m(iCol) = mean(vec);
end;
