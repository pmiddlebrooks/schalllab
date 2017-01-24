function R = parse_JuSu(fn, varnames, trialSelection)

fid = fopen(fn, 'r');
varline = fgetl(fid);
fclose(fid);
varnames0 = regexp(varline, '\t', 'split');
for iVar = 1:length(varnames0),
    f = strfind(varnames0{iVar}, '[');
    if ~isempty(f),
        f2 = strfind(varnames0{iVar}, ']');
        varnames0{iVar}(f:f2) = [];
    end;
    estr = [varnames0{iVar} ' = ' num2str(iVar) ';'];
    eval(estr);
end;

D = dlmread(fn, '\t', 1, 0);

R.varnames = varnames;
indeps = [];
for n = 1:length(R.varnames),
    eval(['indeps = [indeps ' R.varnames{n} '];']);
end;

for iTS = 1:length(trialSelection),
    commandstr = ['vec = D(:, ' trialSelection{iTS}{1} ');'];
    eval(commandstr);
    f = find(vec < trialSelection{iTS}{2} | vec > trialSelection{iTS}{3});
    D(f, :) = [];
end;

levels = [];
levelvals = {};
for n = 1:length(indeps),
    u = unique(D(:, indeps(n)));
    levels = [levels length(u)];
    levelvals{n} = u(:);
end;
R.levels = levels;
R.levelvals = levelvals;

facc = find(D(:, acc) == 1);
dep = D(facc, RT);
indep = D(facc, indeps);
[R.RT, R.idvec, Nvec] = rec_combo(dep, indep);
R.RT = [R.RT Nvec];

dep = D(:, acc);
indep = D(:, indeps);
[R.acc, dum, Nvec] = rec_combo(dep, indep);
R.acc = [R.acc Nvec];
