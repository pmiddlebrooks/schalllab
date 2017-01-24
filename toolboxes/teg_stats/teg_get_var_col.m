function index = teg_get_var_col(varname, varnames)

index = 0;
for ic = 1:length(varnames),
    if strcmp(varname, varnames{ic}),
        index = ic;
        break;
    end;
end;
