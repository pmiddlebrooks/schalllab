function transform_spike_lfp_data(subject)

localDataPath = ['/Users/paulmiddlebrooks/Dropbox/local_data/',lower(subject),'/'];

d = dir(localDataPath);

% for i = 1 : size(d, 1)
for i = 208 : 208
    spikeFlag = 0;
    lfpFlag = 0;
    if regexp(d(i).name, '.*n0.*.mat')
        tic
        disp(d(i).name(1:end-4))
        
        % Check to see if the file has spike or lfp data. If it does, splay out
        % the data, naming each column by its ID
        load(fullfile(localDataPath,d(i).name),'SessionData')
        if isfield(SessionData, 'spikeUnitArray') || isfield(SessionData, 'lfpChannel')
            
            trialData = load(fullfile(localDataPath,d(i).name));
            
            % Transform the spike data
            if isfield(trialData, 'spikeData')
                
                spikeFlag = 1;
                for jUnit = 1 : length(trialData.SessionData.spikeUnitArray)
                    trialData.(trialData.SessionData.spikeUnitArray{jUnit}) = trialData.spikeData(:,jUnit);
                end
                trialData = rmfield(trialData, 'spikeData');
                clear spikeData
            end
            
            
            % Transform the lfp data
            if isfield(trialData, 'lfpData')
                
                lfpFlag = 1;
                
                for jUnit = 1 : length(trialData.SessionData.lfpChannel)
                    unitName = sprintf('lfp%s', num2str(trialData.SessionData.lfpChannel(jUnit), '%02i')); %figure out the channel name
                    trialData.(unitName) = trialData.lfpData(:,jUnit);
                end
                trialData = rmfield(trialData, 'lfpData');
                clear lfpData
            end
        end
        
        if spikeFlag || lfpFlag
            save(fullfile(local_data_path, subject, d(i).name(1:end-4)), '-struct', 'trialData','-v7.3')
        end
        toc
    end
end





