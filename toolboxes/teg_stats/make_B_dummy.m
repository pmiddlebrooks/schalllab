function Bdum = make_B_dummy(BM)

nLevels = length(unique(b));
[X1, factorStarts1, nColsfactor1, labels, cellsets, factors] = teg_create_ANOVA_dummy(BM);

