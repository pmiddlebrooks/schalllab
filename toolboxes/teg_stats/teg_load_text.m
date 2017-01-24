function O = teg_load_text(fn)

% function O = teg_load_text(fn)
% 
%  Assumes first line = variable names.

fid = fopen(fn, 'r');
l0 = fgetl(fid);
fclose(fid);

varnames = regexp(l0, '\t', 'split');
for n = 1:length(varnames),
    if isempty(varnames{n}),
        varnames{n} = ['empty' num2str(n)];
    end;
    eval(['O.vc.' varnames{n} ' = ' num2str(n) ';']);
    eval(['O.varnames{' num2str(n) '} = ''' varnames{n} ''';']);
end;

O.D = dlmread(fn, '\t', 1, 0);
