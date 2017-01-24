% multi-dimensional analysis with SDFs
% make sure you are in proper directory for the alignment you want

%clear;

cd('~/Desktop/zackdata/SDFs/Sac1');

files = dir('S5*');

for i=1:length(files)
    eval(['load ' files(i).name]);
    allmat(i,:,:) = cell2mat(allsdf);
end

% pgm: Get it in the right shape
remat = permute(allmat, [1 3 2]);
remat(isnan(remat)) = 0;
[n t c] = size(remat);

%ignore
longmat = [remat(:,:,1) remat(:,:,2) remat(:,:,3) remat(:,:,4)];

% pgm: 2 X 2 because there are 2 conditions with 2 possible outcomes each
lastmat = reshape(remat, [n,t,2,2]);

% ignore
TT = [0.*ones(1,n) 0.*ones(1,n) 1.*ones(1,n) 1.*ones(1,n)];
RR = [0.*ones(1,n) 1.*ones(1,n) 0.*ones(1,n) 1.*ones(1,n)];

% for r = 1:t;
%      holder = squeeze(remat(:,r,:,:));
% %     newh = [hold(:,1); hold(:,2); hold(:,3); hold(:,4)];
% %     p(r, :) = anovan(newh, {TT RR}, 'display', 'off');%, 'model', 'interaction');
%     pS0_S1(r) = ttest2(holder(:,1), holder(:,2));
%     pI0_I1(r) = ttest2(holder(:,3), holder(:,4));
% end


% pgm: Zach says don't use more than 1/3 of total neurons
np = 6; % # of components
% pgm: Original verision doesn't outpout covs, zach added it
% pgm: black : self-selected
% pgm: red : computer selected
% pgm: solide: rule1
% pgm: dashed: rule2
[W, covs] = dpca(lastmat, np, [], []);
for p=1:np
    Z = W(:,np+1-p)'*longmat; % pgm: np+1-p is to order the plots 
    figure(p)
    plot(1:1401, Z(1:1401), 'k-', ...
       1:1401, Z(1402:2802), 'k-.', ...
       1:1401, Z(2803:4203), 'r-', ...
       1:1401, Z(4204:5604), 'r-.')
end

% Calculating variance of Y (longmat) and marginals
means = mean(longmat')';
meansMat = repmat(means, [1 5604]);
V = sum(sum((longmat-meansMat).^2))./numel(longmat);
% or
Cov = sum(diag(cov(longmat')));
for n=1:7;
    lilCov(n) = sum(diag(covs{n}));
end

%% Calculating explained variances
% Full Y, each PC
% pgm: x2 : each row is a principal component (conditions)
% pgm: x2: each column is the source of the variance (marginals and
% interactions): e.g. [time rule1 rule2 time-rule1-interaction...]
for k=1:np;
    x1(k) = W(:,k)'*cov(longmat')*W(:,k);
    for j=1:7;
        x2(k,j) = W(:,k)'*covs{j}*W(:,k);
    end
end



% Test out explained variance of PCs
Dall = zeros(np, 78, 78);
for k=1:np;
    D(:,:) = W(:,k)*W(:,k)';
    if k==1;
        Dall(k,:,:) = D(:,:);
    else Dall(k,:,:) = squeeze(Dall(k-1,:,:)) + D(:,:);
    end
    pc_longmat = D*longmat;
    lessCov(k) = sum(diag(cov(longmat')));
end