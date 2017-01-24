function [p, ab, ci] = teg_mediation(x, m, y)

nIts = 100;

x = x - mean(x);
m = m - mean(m);
% y = y - mean(y);

if var(x) == 0 || var(m) == 0,
    ab = [NaN NaN];
    p = NaN;
    ci = [NaN NaN];
    return;
end;

[a0, b0] = inner_teg_mediation(x, m, y);
ab0 = a0 * b0;
ab = [a0 b0];

if isnan(b0),
    p = NaN;
    ci = [NaN NaN];
    return;
end;

abv = [];
for iIt = 1:nIts,
    mcf = 1 + floor(rand(size(x)) * length(x));
    xmc = x(mcf);
    mmc = m(mcf);
    ymc = y(mcf);
    [a1, b1] = inner_teg_mediation(xmc, mmc, ymc);
    abv = [abv; a1 * b1];
end;

sgnabv = sign(ab0);
f = find(sgnabv * abv < 0);
p = length(f) / nIts;

abv = sort(abv);
lim1 = floor(0.05 * length(abv));
lim2 = ceil(0.95 * length(abv));
ci = abv([lim1 lim2]);

function [a0, b0] = inner_teg_mediation(x, m, y)
X = x;
Y = m;
a0 = inv(X' * X) * X' * Y;

X = [x m];
Y = y;
[p, t, df, b] = teg_BLR(X, y);
b0 = b(3);

% X = [m];
% Y = y;
% [p, t, df, b] = teg_BLR(X, y);
% b0 = b(2);
