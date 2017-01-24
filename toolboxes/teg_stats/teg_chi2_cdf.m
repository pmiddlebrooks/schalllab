function p = teg_chi2_cdf(x, k)
p = 1 - gammainc(x / 2, k / 2);
