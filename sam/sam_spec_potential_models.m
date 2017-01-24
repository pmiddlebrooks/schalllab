function models = sam_spec_potential_models(SAM)
% SAM_SPEC_POTENTIAL_MODELS Specifies all possible models, given the
% free parameters
%  
% DESCRIPTION 
% <Describe more extensively what this function does> 
%  
% SYNTAX 
% models = SAM_SPEC_POTENTIAL_MODELS(features,nStm,nRsp,nCnd); 
%  
% EXAMPLES
% features = false(3,9,2);  % 3 factors, 9 param. categories, 2 classes
% features(3,1:5,1) = 1;    % 5 go parameters vary across conditions:
%                           %  z0, zc, v, ve, and t0
% nStm = [2 1];             % 2 go-stimuli, 1 stop-stimulus
% nRsp = [2 1];             % go is choice response, stop is simple response
% nCnd = 3;                 % go and stop trials are presented in 3 conditions
%
% models = SAM_SPEC_POTENTIAL_MODELS(features,nStm,nRsp,nCnd);
% 
% REFERENCES 
%  
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Tue 21 Jan 2014 11:47:43 CST by bram 
% $Modified: Tue 21 Jan 2014 11:47:43 CST by bram 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1. PROCESS INPUTS AND DEFINE VARIABLES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% 1.1. Process inputs
% =========================================================================

features = SAM.model.features;

% 1.2. Dynamic variables
% =========================================================================

% Set up models structure
% ------------------------------------------------------------------------- 
models = struct('i',[], ...
                'features',false(size(features)), ...
                'parents',[], ...
                'XSpec',[]);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2. SPECIFY POTENTIAL MODELS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% Number of features of the model (i.e. the number of ways that parameters
% can vary across experimental factors, see Donkin, Chris, Scott Brown, and 
% Andrew Heathcote. ?Drawing Conclusions from Choice Response Time Models: 
% A Tutorial Using the Linear Ballistic Accumulator.? Journal of 
% Mathematical Psychology 55, no. 2 (2011): 140?151).

% Number of features
nFeatures   = sum(features(:));

% Feature indices
iFeature    = find(features);

% Number of possible models (feature combinations)
nModel = nan(1,nFeatures); 
for iF = 1:nFeatures
   nModel(iF) = size(nchoosek(1:nFeatures,iF),1);
end
nModel = sum(nModel) + 1;    % N.B. We add 1 to account for the model in which none of the parameters vary across experimental factors

% Update size of models structure
models = repmat(models,nModel,1);

% Fill in model indices
iCell = num2cell(1:nModel);
[models(1:nModel).i] = deal(iCell{:});

% Specify each model
iModel = 1;
for iF = 1:nFeatures
   combis = nchoosek(1:nFeatures,iF);
   for iCombi = 1:size(combis,1)
      iModel = iModel + 1;
      models(iModel).features(iFeature(combis(iCombi,:))) = 1;
   end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 3. IDENTIFY HIERARCHICAL RELATIONSHIPS AMONG POTENTIAL MODELS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

for iModel = 1:nModel
  parents = false(1,nModel);
  for jModel = 1:nModel
    if jModel ~= iModel
      parents(1,jModel) = all(ismember(find(models(jModel).features(:)),find(models(iModel).features(:))));
    end
  end  
  models(iModel).parents = find(parents);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 4. SPECIFY PARAMETER DETAILS (NUMBER, STATES, INDICES, NAMES)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

for iModel = 1:nModel
  models(iModel).XSpec = sam_spec_x(SAM,models(iModel).features);
end