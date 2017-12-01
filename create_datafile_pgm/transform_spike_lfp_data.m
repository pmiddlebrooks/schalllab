function transform_spike_lfp_data(subject)

localDataPath = ['~/Dropbox/local_data/',lower(subject),'/'];
tebaDataPath = '/Volumes/SchallLab/data/';

switch lower(subject)
    case 'joule'
        fileName = [subject, '.mat'];
        tebaDataPath = [tebaDataPath, 'Joule/'];
    case 'broca'
        fileName = [subject, '.mat'];
        tebaDataPath = [tebaDataPath, 'Broca/'];
    case 'xena'
        fileName = [subject, '.mat'];
        tebaDataPath = [tebaDataPath, 'Xena/Plexon/'];
    otherwise
        fprintf('%s is not a valid subject ID, try again?/n', subject)
        return
end

% d = dir(localDataPath);
d = dir(tebaDataPath);






% % //////////////////////////////////////////////////////////////////////
% %   Convert spikeData and lfpData into single variables for each spike and
% %   lfp channels
% % //////////////////////////////////////////////////////////////////////
%
% for i = 1 : size(d, 1)
%     spikeFlag = 0;
%     lfpFlag = 0;
%     i
%
%     if regexp(d(i).name, '.*n0.*.mat')
%         tic
%
%         % Check to see if the file has spike or lfp data in the old format. If it does, splay out
%         % the data, naming each column by its ID
%         v = who(matfile(fullfile(tebaDataPath,d(i).name)), 'spikeData', 'lfpData');
%
%         if ~isempty(v)
%             disp(d(i).name(1:end-4))
%
%
%             trialData = load(fullfile(tebaDataPath,d(i).name));
%
%             % Transform the spike data
%             if isfield(trialData, 'spikeData')
%
%                 spikeFlag = 1;
%                 for jUnit = 1 : length(trialData.SessionData.spikeUnitArray)
%                     trialData.(trialData.SessionData.spikeUnitArray{jUnit}) = trialData.spikeData(:,jUnit);
%                 end
%                 trialData = rmfield(trialData, 'spikeData');
%                 clear spikeData
%             end
%
%
%             % Transform the lfp data
%             if isfield(trialData, 'lfpData')
%
%                 lfpFlag = 1;
%
%                 for jUnit = 1 : length(trialData.SessionData.lfpChannel)
%                     unitName = sprintf('lfp%s', num2str(trialData.SessionData.lfpChannel(jUnit), '%02i')); %figure out the channel name
%                     trialData.(unitName) = trialData.lfpData(:,jUnit);
%                 end
%                 trialData = rmfield(trialData, 'lfpData');
%                 clear lfpData
%             end
%
%             if spikeFlag || lfpFlag
%                 %                 save(fullfile(local_data_path, subject, d(i).name(1:end-4)), '-struct', 'trialData','-v7.3')
%                 % Save to teba also
%                 save(fullfile(tebaDataPath, d(i).name(1:end-4)), '-struct', 'trialData','-v7.3')
%             end
%             clear trialData
%         end
%
%         toc
%     end
% end







% //////////////////////////////////////////////////////////////////////
%   Create separate files with the lfp data to save memory/load times in
%   python for later when analyzing spiking data
% //////////////////////////////////////////////////////////////////////
for i = 1 : size(d, 1)
    i
    
    if ~isempty(regexp(d(i).name, '.*n0.*.mat')) && isempty(regexp(d(i).name, '.*lfp.*.mat')) && isempty(regexp(d(i).name, '.*legacy.*.mat')) && ~strncmp(d(i).name, '._', 2)
        tic
        
        % Check to see if the file has lfp data. If it does, save an extra
        % file with only the lfp data, and save the old file without the
        % lfp data.
        %         var = who(matfile(fullfile(localDataPath,d(i).name)));
        var = who(matfile(fullfile(tebaDataPath,d(i).name)));
        lfpInd = find(strncmp(var, 'lfp', 3));
        
        % if the file has lfp data
        if ~isempty(lfpInd)
            disp(d(i).name(1:end-4))
            
            %             trialData = load(fullfile(local_data_path, subject ,d(i).name));
            trialData = load(fullfile(tebaDataPath ,d(i).name));
            
            lfpData = struct();
            for j = 1: length(lfpInd)
                jInd = lfpInd(j);
                lfpData.(var{jInd}) = trialData.(var{jInd});
                trialData = rmfield(trialData, var{jInd});
            end
            %         save(fullfile(local_data_path, subject, d(i).name(1:end-4)), '-struct', 'trialData','-v7.3')
            %         save(fullfile(local_data_path, subject, [d(i).name(1:end-4),'_lfp']), '-struct', 'lfpData','-v7.3')
            save(fullfile(local_data_path, subject, d(i).name(1:end-4)), '-struct', 'trialData')
            save(fullfile(local_data_path, subject, [d(i).name(1:end-4),'_lfp']), '-struct', 'lfpData')
            % Save to teba also
            %         save(fullfile(tebaDataPath, d(i).name(1:end-4)), '-struct', 'trialData','-v7.3')
            %         save(fullfile(tebaDataPath, [d(i).name(1:end-4),'_lfp']), '-struct', 'lfpData','-v7.3')
            save(fullfile(tebaDataPath, d(i).name(1:end-4)), '-struct', 'trialData')
            save(fullfile(tebaDataPath, [d(i).name(1:end-4),'_lfp']), '-struct', 'lfpData')
            clear trialData lfpData
        end
        
        toc
        
    end
    
end




% //////////////////////////////////////////////////////////////////////
%   Transform tables to variables in save files
% //////////////////////////////////////////////////////////////////////

% subject = 'joule'
% % dataDir = '/Volumes/SchallLab/data/Joule';
% tebaDataDir = ['/Volumes/SchallLab/data/', subject];
%
% d = dir(fullfile(local_data_path, subject));
% for i = 1 : size(d, 1)
% % for i = 12 : 30
% % for i = 209 : 209
%    if regexp(d(i).name, 'jp.*mat')
%         tic
%           disp(i)
%         disp(d(i).name(1:end-4))
%
%         clear trialOutcome
%
%         load(fullfile(local_data_path, subject ,d(i).name), 'trialOutcome');
%                 if ~exist('trialOutcome')
%
%         load(fullfile(local_data_path, subject ,d(i).name));
% %         load(fullfile(dataDir, d(i).name))
%         trialData = table2struct(trialData, 'ToScalar',true);
%         trialData.SessionData = SessionData;
%
%         save(fullfile(tebaDataDir, d(i).name), '-struct', 'trialData','-v7.3')
%         save(fullfile(local_data_path, subject, d(i).name), '-struct', 'trialData','-v7.3')
%         clear trialData SessionData
%                 end
%         disp(toc)
%     end
% end

% %%  Transform tables to variables in save files
% subject = 'broca'
% % dataDir = '/Volumes/SchallLab/data/Joule';
% tebaDataDir = ['/Volumes/SchallLab/data/', subject];
%
% d = dir(fullfile(local_data_path, subject));
% for i = 1 : size(d, 1)
%     if regexp(d(i).name, 'bp.*mat')
%         tic
%         disp(i)
%         disp(d(i).name(1:end-4))
%
%         clear trialOutcome
%
%         load(fullfile(local_data_path, subject ,d(i).name), 'trialOutcome');
%         if ~exist('trialOutcome')
%
%             load(fullfile(local_data_path, subject ,d(i).name));
%             %         load(fullfile(dataDir, d(i).name))
%             trialData = table2struct(trialData, 'ToScalar',true);
%             trialData.SessionData = SessionData;
%
%             save(fullfile(tebaDataDir, d(i).name), '-struct', 'trialData','-v7.3')
%             save(fullfile(local_data_path, subject, d(i).name), '-struct', 'trialData','-v7.3')
%             clear trialData SessionData
%         end
%         disp(toc)
%     end
% end
%
%
% %%  Transform tables to variables in save files
% subject = 'xena'
% % dataDir = '/Volumes/SchallLab/data/Joule';
% tebaDataDir = ['/Volumes/SchallLab/data/', subject];
%
% d = dir(fullfile(local_data_path, subject));
% for i = 1 : size(d, 1)
%     if regexp(d(i).name, 'xp.*mat')
%         tic
%         disp(i)
%         disp(d(i).name(1:end-4))
%
%         clear trialOutcome
%
%         load(fullfile(local_data_path, subject ,d(i).name), 'trialOutcome');
%         if ~exist('trialOutcome')
%
%             load(fullfile(local_data_path, subject ,d(i).name));
%             %         load(fullfile(dataDir, d(i).name))
%             trialData = table2struct(trialData, 'ToScalar',true);
%             trialData.SessionData = SessionData;
%
%             save(fullfile(tebaDataDir, d(i).name), '-struct', 'trialData','-v7.3')
%             save(fullfile(local_data_path, subject, d(i).name), '-struct', 'trialData','-v7.3')
%             clear trialData SessionData
%         end
%         disp(toc)
%     end
% end



