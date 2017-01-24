function [p, ab, ci] = test_teg_mediation

N = 250;
sc1 = 0.01;
sc2 = 0.01;

x = floor(rand(N, 1) * 2);
m = 4 * x + sc1 * randn(N, 1);
% y = 0.75 * m + sc2 * randn(N, 1);
% y = x + sc2 * randn(N, 1);
y = rand(N, 1);

% figure;
% subplot(3, 1, 1);
% scatter(x, m);
% subplot(3, 1, 2);
% scatter(m, y);
% subplot(3, 1, 3);
% scatter(x, y);

[p, ab, ci] = teg_mediation(x, m, y)
