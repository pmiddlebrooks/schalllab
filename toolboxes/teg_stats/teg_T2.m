function [T2, F, signi, T2Mat] = teg_T2(M)

% function [T2, F, p, T2Mat] = teg_T2(M)
%
% M: variables in columns

C = cov(M);
m = mean(M);
[n, p] = size(M);
T2 = n * m * inv(C) * m(:);
F = T2 * (n - p) / (p * (n - 1));
signi = teg_fsig(F, p, n - p);

M0 = M;
T2Mat = [];
for iObs = 1:n,
    M = M0;
    M(iObs, :) = [];
    o = M0(iObs, :);
    C = cov(M);
    p = size(M, 2);
    T2 = n * (o - m) * inv(C) * (o(:) - m(:));
    T2Mat = [T2Mat; T2];
end;
