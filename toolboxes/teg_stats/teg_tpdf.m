function p = teg_tpdf(t, df)

x = (t + sqrt(t.^2 + df)) / (2 * sqrt(t.^2 + df));
z = df / 2;
w = df / 2;
tcdf00 = betainc(x, z, w);
p = 1 - tcdf00;
if p > 0.5,
    p = 1 - p;
end;
p = p * 2;
