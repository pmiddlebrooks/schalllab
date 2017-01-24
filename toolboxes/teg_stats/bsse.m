function se = bsse(M)

% Bootstrapped SE per column

nIts = 5;
bsM = [];
for iIt = 1:nIts,
    randsel = 1 + floor(size(M, 1) * rand(size(M, 1), 1));
    bsM0 = M(randsel, :);
    m0 = mean(bsM0);
    bsM = [bsM; bsM0];
end;

m = mean(bsM);
se = sqrt(var(bsM)); % SD of the mean
