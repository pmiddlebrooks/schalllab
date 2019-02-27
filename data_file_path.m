function [tebaDataFile, localDataPath, localDataFile] = data_file_path(subjectID, sessionID, monkeyOrHuman)

if nargin < 3
if ismember(lower(subjectID), {'joule', 'broca', 'xena', 'chase', 'hoagie', 'norm', 'andy','shuffles','nebby'})
   monkeyOrHuman = 'monkey';
else
      monkeyOrHuman = 'human';
end
end


homeDataPath = '/Users/paulmiddlebrooks/schalllab/local_data/';
homeDataPath = '/Users/paulmiddlebrooks/Dropbox/local_data/';
tebaDataPath = '/Volumes/SchallLab/data/';
if isdir(tebaDataPath)
   location = 'work';
elseif isdir(homeDataPath)
   location = 'home';
else
   disp('If you''re at work you may need to connect to teba')
   return
end

switch monkeyOrHuman
   case 'monkey'
      
%       localDataPath = ['/Users/paulmiddlebrooks/schalllab/local_data/',lower(subjectID),'/'];
      localDataPath = ['/Users/paulmiddlebrooks/Dropbox/local_data/',lower(subjectID),'/'];
      switch location
         case 'work'
            
            switch lower(subjectID)
               case 'joule'
                  fileName = [sessionID, '.mat'];
                  tebaDataFile = [tebaDataPath, 'Joule/', fileName];
               case 'broca'
                  fileName = [sessionID, '.mat'];
                  tebaDataFile = [tebaDataPath, 'Broca/', fileName];
               case 'xena'
                  fileName = [sessionID, '.mat'];
                  tebaDataFile = [tebaDataPath, 'Xena/Plexon/', fileName];
               case 'andy'
                  fileName = [sessionID, '.mat'];
                  tebaDataFile = [tebaDataPath, 'andy/andyfef/PDP', fileName];
               case 'chase'
                  fileName = [sessionID, '.mat'];
                  tebaDataFile = [tebaDataPath, 'chase/chafef/pdp', fileName];
               case 'nebby'
                  fileName = [sessionID, '.mat'];
                  tebaDataFile = [];
               case 'shuffles'
                  fileName = [sessionID, '.mat'];
                  tebaDataFile = [];
               otherwise
                  fprintf('%s is not a valid subject ID, try again?/n', subjectID)
                  return
            end
            
         case 'home'
            tebaDataFile = [];
            switch lower(subjectID)
               case 'joule'
                  fileName = [sessionID, '.mat'];
               case 'broca'
                  fileName = [sessionID, '.mat'];
               case 'xena'
                  fileName = [sessionID, '.mat'];
               case 'nebby'
                  fileName = [sessionID, '.mat'];
               case 'shuffles'
                  fileName = [sessionID, '.mat'];
               case 'andy'
                  fileName = [sessionID, '.mat'];
               case 'chase'
                  fileName = [sessionID, '.mat'];
               otherwise
                  fprintf('%s is not a valid subject ID, try again?\n', subjectID)
                  return
            end
            
      end
      
      
      
      
      
      
      
   case 'human'
      
      humanDataPath = '/Volumes/middlepg/HumanData/ChoiceStopTask/';
      
      localDataPath = '/Users/paulmiddlebrooks/schalllab/local_data/human/';
      switch location
         case 'work'
            
            fileName = [sessionID, '.mat'];
            tebaDataFile = [humanDataPath, fileName];
            
         case 'home'
            tebaDataFile = [];
            fileName = [sessionID, '.mat'];
      end
      
end
localDataFile = fullfile(localDataPath, fileName);


end