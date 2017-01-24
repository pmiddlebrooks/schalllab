function [modelvec, dirpar] = teg_trendy(vec)

vec = vec(:);
dirpar = 0;

f = find(isnan(vec));
if ~isempty(f),
    ff = find(~isnan(vec));
    vec(f) = mean(vec(ff));
end;

modelvec = vec;
dirpar = mean(diff(vec));
return;

% newvec = zeros(size(vec));
% for iVec = 1:length(vec),
%     a = max(1, iVec - 1);
%     b = min(length(vec), iVec + 1);
%     newvec(iVec) = mean(vec(a:b));
% end;
% modelvec = newvec;

t = 1:length(vec);
t = t - 1;
t = t ./ max(t);
t = t(:);
X = [];
f = 0.5;
x0 = ones(size(t));
X = [X x0(:)];
% x0 = t - mean(t);
% X = [X x0(:)];
x0 = cos(2 * pi * f * t);
X = [X x0(:)];
b = inv(X'*X)*X'*vec;
modelvec = X * b;
dirpar = mean(b(2:end));

% t = 1:length(vec);
% t = t - 1;
% t = t ./ max(t);
% t = t * 2 * pi;
% X = [];
% for iComp = 0:1,
%     x0 = t .^ iComp;
%     X = [X x0(:)];
% end;
% b = inv(X'*X)*X'*vec;
% modelvec = X * b;
% dirpar = b(2);
