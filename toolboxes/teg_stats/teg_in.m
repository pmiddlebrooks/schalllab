function in0 = teg_in(v1, v2)

v1 = v1(:)';
v2 = v2(:);

M1 = ones(length(v2), 1) * v1;
M2 = v2 * ones(1, length(v1));

Comp0 = M1 == M2;

in0 = length(find(Comp0(:))) > 0;
