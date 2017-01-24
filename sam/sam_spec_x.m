function XSpec = sam_spec_x(SAM,features)
% SAM_SPEC_X <Synopsis of what this function does> 
%  
% DESCRIPTION 
% <Describe more extensively what this function does> 
%  
% SYNTAX 
% SAM_SPEC_X; 
%  
% EXAMPLES 
%  
%  
% REFERENCES 
%  
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Wed 29 Jan 2014 13:12:39 CST by bram 
% $Modified: Wed 29 Jan 2014 13:12:39 CST by bram 
 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1. PROCESS INPUTS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 
% % Static variables
% ========================================================================= 

nStruct     = struct('n',[], ...
                     'nCat',[], ...
                     'nCatClass',[]);
freeStruct  = struct('free',[], ...
                     'freeCat',[], ...
                     'freeCatClass',[]);
iStruct     = struct('iCat',[], ...
                     'iCatClass',[]);
nameStruct  = struct('name',[], ...
                     'nameCat',[], ...
                     'nameCatClass',[]);               
classStruct = struct('go',[],...
                     'stop',[],...
                     'all',[]);
                   XSpec       = struct('nCombi',[], ...
                     'n',nStruct, ...
                     'free',classStruct, ...
                     'i',classStruct, ...
                     'name',classStruct);

% % Dynamic variables
% ========================================================================= 

nClass        = size(features,3);

nStm          = SAM.expt.nStm;
nRsp          = SAM.expt.nRsp;
nCnd          = SAM.expt.nCnd;

maxNStm       = max(cell2mat(nStm(:)),[],1);
maxNRsp       = max(cell2mat(nRsp(:)),[],1);

taskFactors   = [maxNStm;maxNRsp;nCnd,nCnd];

included      = SAM.model.XCat.included;

classSpecific = SAM.model.XCat.classSpecific;

nXCat         = SAM.model.XCat.n;

XCatName      = SAM.model.XCat.name;

className     = SAM.model.general.classNames;

iVe           = SAM.model.XCat.i.iVe;

iWliw         = SAM.model.XCat.i.iWliw;

iWffiw        = SAM.model.XCat.i.iWffiw;

iWlib         = SAM.model.XCat.i.iWlib;

iScale        = SAM.model.XCat.scale.iX;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2. DETERMINE ALL FACTORIAL COMBINATIONS OF TASK FACTORS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

nCombi = cell(1,nClass);

for iClass = 1:nClass
  nCombi{iClass} = diag(taskFactors(:,iClass)) * any(features(:,:,iClass),2);
  nCombi{iClass}(nCombi{iClass} == 0) = 1;
  nCombi{iClass} = prod(nCombi{iClass});
end

% Put variables in output structure
XSpec.nCombi  = nCombi;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 3. DETERMINE THE NUMBER OF PARAMETERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Number per task factor, per parameter category, per accumulator class
nCat = zeros(size(features(:,:,1:nClass)));
for iClass = 1:nClass
  nCat(:,:,iClass) = diag(taskFactors(:,iClass)) * features(:,:,iClass);
end
nCat(nCat == 0) = 1;

% Per parameter category, per accumulator class
nCatClass = prod(nCat,1);

% Per parameter category
nCat = sum(nCatClass,3);

% Correct for excluded parameters
nCat(~included) = 1;
nCatClass(:,~included,:) = 1;

% Correct for parameters that are not class-specific
nCat(~classSpecific) = 1;
nCatClass(:,~classSpecific) = 1;

% Correct Ve, Wliw, and Wffiw for classes that have only one accumulator
if included(iVe)
  nCatClass(:,iVe,maxNRsp(1:nClass) <= 1) = 0;
  nCat(iVe) = sum(nCatClass(:,iVe,:),3);
end

if included(iWliw)
  nCatClass(:,iWliw,maxNRsp(1:nClass) <= 1) = 0;
  nCat(iWliw) = sum(nCatClass(:,iWliw,:),3);
end

if included(iWffiw)
  nCatClass(:,iWffiw,maxNRsp(1:nClass) <= 1) = 0;
  nCat(iWffiw) = sum(nCatClass(:,iWffiw,:),3);
end

% % Constrain lateral inhibition between classes to effects of STOP onto GO
% % only, not vice versa
% if ~included(iWlib)
%   nCatClass(:,iWlib,1) = 0;
%   nCat(iWlib) = sum(nCatClass(:,iWlib,:),3);
% end

% Squeeze out redundant dimensions, if any
nCatClass = reshape(nCatClass,nXCat,nClass)';

% Put variables in output structure
XSpec.n.n                = sum(nCat);
XSpec.n.nCat             = nCat;
XSpec.n.nCatClass        = nCatClass;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 4. DETERMINE THE STATE OF PARAMETERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Identify free parameters, per category
freeCat = arrayfun(@(a,b) repmat(a,1,b),included,nCat,'Uni',0);

% Correct for the scaling parameter
freeCat{iScale} = false;

freeCatClass = cell(nClass,nXCat);

for iXCat = 1:nXCat
  
  if ~included(iXCat)
    % If: 
    % - parameter category not included,
    
    freeCatClass(:,iXCat) = repmat(freeCat(iXCat),nClass,1);
    
  elseif included(iXCat) && ~classSpecific(iXCat)
    % If:
    % - parameter category included, 
    % - not class-specific
    
    freeCatClass(:,iXCat) = repmat(freeCat(iXCat),nClass,1);
  elseif included(iXCat) && classSpecific(iXCat)
    % If:
    % - parameter category included, 
    % - class-specific
    
    freeCatClass(:,iXCat) = mat2cell(freeCat{iXCat}(:),nCatClass(:,iXCat),1);
  end
end

% Ensure row vectors
freeCatClass = cellfun(@(in1) in1(:)',freeCatClass,'Uni',0);

% Set lateral inhibition from GO to STOP to 0
freeCatClass{1,iWlib} = false(size(freeCatClass{1,iWlib}));
iWlibGO = find(~freeCatClass{1,iWlib});
freeCat{iWlib}(iWlibGO) = false;

% Put variables in output structure
% =========================================================================

% GO parameters
% -------------------------------------------------------------------------
XSpec.free.go                 = freeStruct;
XSpec.free.go.free            = cell2mat(freeCatClass(1,:));
XSpec.free.go.freeCat         = freeCatClass(1,:);
XSpec.free.go.freeCatClass    = freeCatClass(1,:);

% STOP parameters
% -------------------------------------------------------------------------
stopFreeCatClass              = freeCatClass;
stopFreeCatClass(1,:)         = cellfun(@(in1) assignit(in1,1:numel(in1),false),stopFreeCatClass(1,:),'Uni',0);

stopFreeCat                   = cell(1,nXCat);

for iXCat = 1:nXCat
  
  if ~included(iXCat)
    % If: 
    % - parameter category not included,
    stopFreeCat{iXCat}        = stopFreeCatClass{2,iXCat};
   
  elseif included(iXCat) && ~classSpecific(iXCat)
    % If:
    % - parameter category included, 
    % - not class-specific
    stopFreeCat{iXCat}        = stopFreeCatClass{2,iXCat};
  elseif included(iXCat) && classSpecific(iXCat)
    % If:
    % - parameter category included, 
    % - not class-specific
    stopFreeCat{iXCat}        = [stopFreeCatClass{:,iXCat}];
  end
end


XSpec.free.stop               = freeStruct;
XSpec.free.stop.free          = [stopFreeCat{:}];
XSpec.free.stop.freeCat       = stopFreeCat;
XSpec.free.stop.freeCatClass  = stopFreeCatClass;

% ALL parameters
% -------------------------------------------------------------------------
XSpec.free.all                = freeStruct;
XSpec.free.all.free           = [freeCat{:}];
XSpec.free.all.freeCat        = freeCat;
XSpec.free.all.freeCatClass   = freeCatClass;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 5. DETERMINE THE INDICES OF PARAMETERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% GO parameters
% -------------------------------------------------------------------------

% First and last index per parameter category
i1    = [1,cumsum(nCatClass(1,1:end-1))+1];
iend  = cumsum(nCatClass(1,:));

% Indices per parameter category
iCat = arrayfun(@(a,b) a:b,i1,iend,'Uni',0);

% iCat and iCatClass are identical. Ensure row vectors
iCatClass = cellfun(@(in1) in1(:)',iCat,'Uni',0);

% Put variables in output structure
XSpec.i.go                = iStruct;
XSpec.i.go.iCat           = iCat;
XSpec.i.go.iCatClass      = iCatClass;

clear i1 iend iCat iCatClass

% STOP parameters
% -------------------------------------------------------------------------

% First and last index per parameter category
i1    = [1,cumsum(nCat(1:end-1))+1];
iend  = cumsum(nCat);

% Indices per parameter category
iCat = arrayfun(@(a,b) a:b,i1,iend,'Uni',0);

iCatClass = cell(nClass,nXCat);

for iXCat = 1:nXCat
  
  if ~included(iXCat)
    % If: 
    % - parameter category not included,
    iCatClass(:,iXCat) = repmat(iCat(iXCat),nClass,1);
  elseif included(iXCat) && ~classSpecific(iXCat)
    % If:
    % - parameter category included, 
    % - not class-specific
    iCatClass(:,iXCat) = repmat(iCat(iXCat),nClass,1);
  elseif included(iXCat) && classSpecific(iXCat)
    % If:
    % - parameter category included, 
    % - not class-specific
    iCatClass(:,iXCat) = mat2cell(iCat{iXCat},1,nCatClass(:,iXCat))';
  end
end


% Ensure row vectors
iCatClass = cellfun(@(in1) in1(:)',iCatClass,'Uni',0);

% Put variables in output structure
XSpec.i.stop              = iStruct;
XSpec.i.stop.iCat         = iCat;
XSpec.i.stop.iCatClass    = iCatClass;

% ALL parameters
% -------------------------------------------------------------------------

% First and last index per parameter category
i1    = [1,cumsum(nCat(1:end-1))+1];
iend  = cumsum(nCat);

% Indices per parameter category
iCat = arrayfun(@(a,b) a:b,i1,iend,'Uni',0);

iCatClass = cell(nClass,nXCat);

for iXCat = 1:nXCat
  
  if ~included(iXCat)
    % If: 
    % - parameter category not included,
    iCatClass(:,iXCat) = repmat(iCat(iXCat),nClass,1);
  elseif included(iXCat) && ~classSpecific(iXCat)
    % If:
    % - parameter category included, 
    % - not class-specific
    iCatClass(:,iXCat) = repmat(iCat(iXCat),nClass,1);
  elseif included(iXCat) && classSpecific(iXCat)
    % If:
    % - parameter category included, 
    % - not class-specific
    iCatClass(:,iXCat) = mat2cell(iCat{iXCat},1,nCatClass(:,iXCat))';
  end
end

% Ensure row vectors
iCatClass = cellfun(@(in1) in1(:)',iCatClass,'Uni',0);

% Put variables in output structure
XSpec.i.all             = iStruct;
XSpec.i.all.iCat        = iCat;
XSpec.i.all.iCatClass   = iCatClass;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 6. DETERMINE THE NAMES OF PARAMETERS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

nameCatClass = cell(nClass,nXCat);

for iXCat = 1:nXCat

  if ~included(iXCat)
    % If: 
    % - parameter category not included,
    nameCatClass(:,iXCat) = repmat({sprintf('%s',XCatName{iXCat})},nClass,1);
  elseif included(iXCat) && ~classSpecific(iXCat)
    % If:
    % - parameter category included, 
    % - not class-specific
    nameCatClass(:,iXCat) = repmat({sprintf('%s',XCatName{iXCat})},nClass,1);
    
  elseif included(iXCat) && classSpecific(iXCat) && nCat(iXCat) == 1 && ~all(nCatClass(:,iXCat))
    fun = @(a) sprintf('%s_{%s}',XCatName{iXCat},a);
    nameCatClass(:,iXCat) = cellfun(fun,className,'Uni',0);
    nameCatClass(~nCatClass(:,iXCat),iXCat) = {''};
  elseif included(iXCat) && classSpecific(iXCat) && nCat(iXCat) == nClass
    fun = @(a) sprintf('%s_{%s}',XCatName{iXCat},a);
    nameCatClass(:,iXCat) = cellfun(fun,className,'Uni',0);
  elseif included(iXCat) && classSpecific(iXCat) && nCat(iXCat) > nClass
    for iClass = 1:nClass
      % Identify how parameter category varies across task factors
      signature = logical(features(:,iXCat,iClass)); 

      % Temporary variables to keep fun readable
      thisCatName = XCatName{iXCat};  % Parameter category name
      thisClassName = className{iClass};   % Class name

      if isequal(signature,[0 0 0]')
        if nCatClass(iClass,iXCat) == 0
          nameCatClass{iClass,iXCat} = '';
        else
          fun = @(a) sprintf('%s_{%s}',thisCatName,a);
          nameCatClass{iClass,iXCat} = cellfun(fun,className(iClass),'Uni',0);
        end
      else

        combi     = fullfact(taskFactors(signature,iClass))';
        nRow      = size(combi,1);
        nCol      = size(combi,2);
        combiCell = mat2cell(combi,nRow,ones(nCol,1));

        if isequal(signature,[1 0 0]')
          fun = @(a) sprintf('%s_{%s,s%d}',thisCatName,thisClassName,a);
        elseif isequal(signature,[0 1 0]')
          fun = @(a) sprintf('%s_{%s,r%d}',thisCatName,thisClassName,a);
        elseif isequal(signature,[0 0 1]')
          fun = @(a) sprintf('%s_{%s,c%d}',thisCatName,thisClassName,a);
        elseif isequal(signature,[1 1 0]')
          fun = @(a) sprintf('%s_{%s,s%d,r%d}',thisCatName,thisClassName,a);
        elseif isequal(signature,[1 0 1]')
          fun = @(a) sprintf('%s_{%s,s%d,c%d}',thisCatName,thisClassName,a);
        elseif isequal(signature,[0 1 1]')
          fun = @(a) sprintf('%s_{%s,r%d,c%d}',thisCatName,thisClassName,a);
        elseif isequal(signature,[1 1 1]')
          fun = @(a) sprintf('%s_{%s,s%d,r%d,c%d}',thisCatName,thisClassName,a);
        end
        nameCatClass{iClass,iXCat} = cellfun(fun,combiCell,'Uni',0);
      end
    end
  end
    
%   if nCat(iXCat) == 1 && all(nCatClass(:,iXCat))
%     nameCatClass(:,iXCat) = repmat({sprintf('%s',XCatName{iXCat})},nClass,1);
%   elseif nCat(iXCat) == 1 && ~all(nCatClass(:,iXCat))
%     fun = @(a) sprintf('%s_{%s}',XCatName{iXCat},a);
%     nameCatClass(:,iXCat) = cellfun(fun,className,'Uni',0);
% %     nameCatClass(:,iXCat) = repmat({sprintf('%s',XCatName{iXCat})},nClass,1);
%     nameCatClass(~nCatClass(:,iXCat),iXCat) = {''};
%   elseif nCat(iXCat) == nClass
%     fun = @(a) sprintf('%s_{%s}',XCatName{iXCat},a);
%     nameCatClass(:,iXCat) = cellfun(fun,className,'Uni',0);
%   elseif nCat(iXCat) > nClass
%     for iClass = 1:nClass
%       % Identify how parameter category varies across task factors
%       signature = logical(features(:,iXCat,iClass)); 
% 
%       % Temporary variables to keep fun readable
%       thisCatName = XCatName{iXCat};  % Parameter category name
%       thisClassName = className{iClass};   % Class name
% 
%       if isequal(signature,[0 0 0]')
%         if nCatClass(iClass,iXCat) == 0
%           nameCatClass{iClass,iXCat} = '';
%         else
%           fun = @(a) sprintf('%s_{%s}',thisCatName,a);
%           nameCatClass{iClass,iXCat} = cellfun(fun,className(iClass),'Uni',0);
%         end
%       else
% 
%         combi     = fullfact(taskFactors(signature,iClass))';
%         nRow      = size(combi,1);
%         nCol      = size(combi,2);
%         combiCell = mat2cell(combi,nRow,ones(nCol,1));
% 
%         if isequal(signature,[1 0 0]')
%           fun = @(a) sprintf('%s_{%s,s%d}',thisCatName,thisClassName,a);
%         elseif isequal(signature,[0 1 0]')
%           fun = @(a) sprintf('%s_{%s,r%d}',thisCatName,thisClassName,a);
%         elseif isequal(signature,[0 0 1]')
%           fun = @(a) sprintf('%s_{%s,c%d}',thisCatName,thisClassName,a);
%         elseif isequal(signature,[1 1 0]')
%           fun = @(a) sprintf('%s_{%s,s%d,r%d}',thisCatName,thisClassName,a);
%         elseif isequal(signature,[1 0 1]')
%           fun = @(a) sprintf('%s_{%s,s%d,c%d}',thisCatName,thisClassName,a);
%         elseif isequal(signature,[0 1 1]')
%           fun = @(a) sprintf('%s_{%s,r%d,c%d}',thisCatName,thisClassName,a);
%         elseif isequal(signature,[1 1 1]')
%           fun = @(a) sprintf('%s_{%s,s%d,r%d,c%d}',thisCatName,thisClassName,a);
%         end
%         nameCatClass{iClass,iXCat} = cellfun(fun,combiCell,'Uni',0);
%       end
%     end
%   end
end

% Correct for inexisting variables (e.g. ve, in classes without
% response alternatives)
if any(nCatClass == 0)
  nameCatClass{nCatClass == 0} = [];
end

nameCat = cell(1,nXCat);

for iXCat = 1:nXCat
  if sum(nCatClass(:,iXCat)) == nCat(iXCat)
    nameCat{iXCat} = getit(nameCatClass(:,iXCat));
    nameCat{iXCat}(cellfun(@isempty,nameCat{iXCat})) = [];
  elseif XSpec.n.nCat(iXCat) == 1
    nameCat{iXCat} = getit(nameCatClass(1,iXCat));
  end
end

% Put variables in output structure
% =========================================================================

% GO parameters
% -------------------------------------------------------------------------
goNames                       = getit(nameCatClass(1,:));
goNames(strcmp('',goNames))   = [];
XSpec.name.go                 = nameStruct;

XSpec.name.go.name            = goNames;
XSpec.name.go.nameCat         = goNames;
XSpec.name.go.nameCatClass    = goNames;

% STOP parameters
% -------------------------------------------------------------------------
XSpec.name.stop               = nameStruct;
XSpec.name.stop.name          = getit(nameCat);
XSpec.name.stop.nameCat       = nameCat;
XSpec.name.stop.nameCatClass  = nameCatClass;

% ALL parameters
% -------------------------------------------------------------------------
XSpec.name.all                = nameStruct;
XSpec.name.all.name           = getit(nameCat);
XSpec.name.all.nameCat        = nameCat;
XSpec.name.all.nameCatClass   = nameCatClass;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7. SUBFUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c, sz] = getit(c)
% Extract data from cell arrays containing cell arrays

 if iscell(c)
     d = size(c,2);
     [c, sz] = cellfun(@getit, c, 'UniformOutput', 0);
     c = cat(2,c{:});
     sz = [sz{1} d];
 else
     c = {c};
     sz = [];
 end
 
function C = assignit(C,i,val)
% Assigns values to indices within a cell
C(i) = val;