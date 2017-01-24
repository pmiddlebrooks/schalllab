function [tebDataFile, localDataPath, localDataFile] = data_file_path(subjectID, sessionID, monkeyOrHuman)

if nargin < 3
if ismember(lower(subjectID), {'joule', 'broca', 'xena', 'chase', 'hoagie', 'norm', 'andy','shuffles','nebby'})
   monkeyOrHuman = 'monkey';
else
      monkeyOrHuman = 'human';
end
end


homeDataPath = '/Users/paulmiddlebrooks/schalllab/local_data/';
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
      
      localDataPath = ['/Users/paulmiddlebrooks/schalllab/local_data/',lower(subjectID),'/'];
      switch location
         case 'work'
            
            switch lower(subjectID)
               case 'joule'
                  fileName = [sessionID, '.mat'];
                  tebDataFile = [tebaDataPath, 'Joule/', fileName];
               case 'broca'
                  fileName = [sessionID, '.mat'];
                  tebDataFile = [tebaDataPath, 'Broca/', fileName];
               case 'xena'
                  fileName = [sessionID, '.mat'];
                  tebDataFile = [tebaDataPath, 'Xena/Plexon/', fileName];
               case 'andy'
                  fileName = [sessionID, '.mat'];
                  tebDataFile = [tebaDataPath, 'andy/andyfef/PDP', fileName];
               case 'chase'
                  fileName = [sessionID, '.mat'];
                  tebDataFile = [tebaDataPath, 'chase/chafef/pdp', fileName];
               case 'nebby'
                  fileName = [sessionID, '.mat'];
                  tebDataFile = [];
               case 'shuffles'
                  fileName = [sessionID, '.mat'];
                  tebDataFile = [];
               otherwise
                  fprintf('%s is not a valid subject ID, try again?/n', subjectID)
                  return
            end
            
         case 'home'
            tebDataFile = [];
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
            tebDataFile = [humanDataPath, fileName];
            
         case 'home'
            tebDataFile = [];
            fileName = [sessionID, '.mat'];
      end
      
end
localDataFile = fullfile(localDataPath, fileName);


end