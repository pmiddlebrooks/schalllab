function [F]=analyze_processedTest(numSubj, plots)
%
%
% plots= 1 or true if plots should be output, zero otherwise
%
% %analyze processed data for GazeDirect
% %Processed Data from ascreaderGazeDirect needs to be in Processed
% Data/subjectnumber_workspace.mat format for all subjects
%
% outputs a struct F of length (number of participants)
% -JI 2012

%*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
%Load all workspaces for all viable subjects
%*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
jj=1;
subjInd=[];

for ii=1:numSubj
    iiNum=num2str(ii);
    if ii<10
        iiNum=['0' iiNum];
    end
    %make sure workspaces are all the same size/don't have these extra
    %fields in them already
    if exist(['ProcessedData/' iiNum '_workspace.mat'])==2
        ll=1;
        tempStruct=importdata(['ProcessedData/' iiNum '_workspace.mat']);
        while ll==1;
            if isfield(tempStruct,'memHits')
                tempStruct=rmfield(tempStruct,'memHits');
            elseif isfield(tempStruct,'memMisses')
                tempStruct=rmfield(tempStruct,'memMisses');
            elseif isfield(tempStruct,'memNoResp')
                tempStruct=rmfield(tempStruct,'memNoResp');
            elseif isfield(tempStruct,'realMemRT')
                tempStruct=rmfield(tempStruct,'realMemRT');
            elseif isfield(tempStruct,'avgMemRT')
                tempStruct=rmfield(tempStruct,'avgMemRT');
            elseif isfield(tempStruct,'realMemAcc')
                tempStruct=rmfield(tempStruct,'realMemAcc');
            elseif isfield(tempStruct,'FirstSacColChg')
                tempStruct=rmfield(tempStruct,'FirstSacColChg');
            elseif isfield(tempStruct,'FirstSacNonObj')
                tempStruct=rmfield(tempStruct,'FirstSacNonObj');    
            elseif isfield(tempStruct,'FirstSacOther')
                tempStruct=rmfield(tempStruct,'FirstSacOther');
            elseif isfield(tempStruct,'FirstSacUnspec')
                tempStruct=rmfield(tempStruct,'FirstSacUnspec');    
            elseif isfield(tempStruct,'percFirstSacColChg')
                tempStruct=rmfield(tempStruct,'percFirstSacColChg');
            elseif isfield(tempStruct,'percFirstSacOther')
                tempStruct=rmfield(tempStruct,'percFirstSacOther');
            elseif isfield(tempStruct,'percFirstSacNonObj')
                tempStruct=rmfield(tempStruct,'percFirstSacNonObj');
             elseif isfield(tempStruct,'percFirstSacUnspec')
                tempStruct=rmfield(tempStruct,'percFirstSacUnspec');
            elseif isfield(tempStruct,'SacColChg')
                tempStruct=rmfield(tempStruct,'SacColChg');
            elseif isfield(tempStruct,'SacColNonObj')
                tempStruct=rmfield(tempStruct,'SacNonObj');
            elseif isfield(tempStruct,'SacOther')
                tempStruct=rmfield(tempStruct,'SacOther');
            elseif isfield(tempStruct,'SacUnspec')
                tempStruct=rmfield(tempStruct,'SacUnspec');    
            elseif isfield(tempStruct,'percSacColChg')
                tempStruct=rmfield(tempStruct,'percSacColChg');
            elseif isfield(tempStruct,'percSacOther')
                tempStruct=rmfield(tempStruct,'percSacOther');
            elseif isfield(tempStruct,'percSacNonObj')
                tempStruct=rmfield(tempStruct,'percSacNonObj'); 
            elseif isfield(tempStruct,'percSacUnspec')
                tempStruct=rmfield(tempStruct,'percSacUnspec');    
            else
                ll=2;
            end
        end
        
        F(jj)=tempStruct;
        jj=jj+1;
        subjInd=[subjInd ii];
    end
end


%*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
%Calculate mean reaction times
%*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
kk=1;
for ii=1:numSubj
    if any(subjInd==ii)
        
        F(1,kk).realMemRT=(F(1,kk).memRT(F(1,kk).memRT>0));
        F(1,kk).avgMemRT=mean(F(1,kk).realMemRT);
        
        %*~*~*~*~*~*~*~*~*~*~*~*~*~*
        %Calculate percent accuracy
        %*~*~*~*~*~*~*~*~*~*~*~*~*~*
        
        F(1,kk).memHits=length(find(F(1,kk).memacc==1));
        F(1,kk).memNoResp=length(find(F(1,kk).memacc==2));
        F(1,kk).memMisses=length(find(F(1,kk).memacc==0));
        F(1,kk).realMemAcc=(F(1,kk).memHits/(F(1,kk).memHits+F(1,kk).memMisses+F(1,kk).memNoResp));
        
        %*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
        %Calculate # & % of first saccades to color change v. other
        %*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~*~*~
        F(1,kk).FirstSacColChg=length(find(F(1,kk).fixcode(1,:)==1));
        F(1,kk).FirstSacOther=length(find(F(1,kk).fixcode(1,:)==2));
        F(1,kk).FirstSacNonObj=length(find(F(1,kk).fixcode(1,:)==3));
        F(1,kk).FirstSacUnspec=length(find(F(1,kk).fixcode(1,:)==4));
        F(1,kk).percFirstSacColChg=(F(1,kk).FirstSacColChg/(F(1,kk).FirstSacColChg+F(1,kk).FirstSacOther+F(1,kk).FirstSacNonObj+F(1,kk).FirstSacUnspec));
        F(1,kk).percFirstSacOther=F(1,kk).FirstSacOther/(F(1,kk).FirstSacColChg+F(1,kk).FirstSacOther+F(1,kk).FirstSacNonObj+F(1,kk).FirstSacUnspec);
        F(1,kk).percFirstSacNonObj=F(1,kk).FirstSacNonObj/(F(1,kk).FirstSacColChg+F(1,kk).FirstSacOther+F(1,kk).FirstSacNonObj+F(1,kk).FirstSacUnspec);
        F(1,kk).percFirstSacUnspec=F(1,kk).FirstSacUnspec/(F(1,kk).FirstSacColChg+F(1,kk).FirstSacOther+F(1,kk).FirstSacNonObj+F(1,kk).FirstSacUnspec);
        %*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
        %Calculate total # and % of saccades to color change v. other
        %*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~*~*~
        F(1,kk).SacColChg=length(find(F(1,kk).fixcode==1));
        F(1,kk).SacOther=length(find(F(1,kk).fixcode==2));
        F(1,kk).SacNonObj=length(find(F(1,kk).fixcode==3));
        F(1,kk).SacUnspec=length(find(F(1,kk).fixcode==4));
        F(1,kk).percSacColChg=(F(1,kk).SacColChg/(F(1,kk).SacColChg+F(1,kk).SacOther+F(1,kk).SacNonObj+F(1,kk).SacUnspec));
        F(1,kk).percSacOther=F(1,kk).SacOther/(F(1,kk).SacColChg+F(1,kk).SacOther+F(1,kk).SacNonObj+F(1,kk).SacUnspec);
        F(1,kk).percSacNonObj=F(1,kk).SacNonObj/(F(1,kk).SacColChg+F(1,kk).SacOther+F(1,kk).SacNonObj+F(1,kk).SacUnspec);
        F(1,kk).percSacUnspec=F(1,kk).SacUnspec/(F(1,kk).SacColChg+F(1,kk).SacOther+F(1,kk).SacNonObj+F(1,kk).SacUnspec);
        kk=kk+1;
    end
end
if plots
    %*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
    %Plot RTs for all subject #s
    %*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~*~*~
    plot([F.avgMemRT],'x')
    title('Reaction Times')
    %*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
    %Plot Percent Accuracy for all subject #s
    %*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~*~*~
    figure
    plot([F.realMemAcc],'x')
    title('Percent Accuracy')
    
    %*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
    %Plot Percent first saccades for all subject #s
    %*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~*~*~
    figure
    plot ([F.percFirstSacColChg],'x')
    hold on
    plot([F.percFirstSacOther],'xr')
    plot([F.percFirstSacNonObj],'xg')
    plot([F.percFirstSacUnspec],'xc')
    title('Percent First Saccades to...')
    legend('Color Change','Other Object','Non-Object','Unspecified')
    hold off
    %*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
    %Plot Percent saccades for all subject #s
    %*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~*~*~
    figure
    plot ([F.percSacColChg],'x')
    hold on
    plot([F.percSacOther],'xr')
    plot([F.percSacNonObj],'xg')
    plot([F.percSacUnspec],'xc')
    title('Percent Saccades to...')
    legend('Color Change','Other Object','Non-Object','Unspecified')
    hold off
       
end

%*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
%Calculate variables of interest 
%*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
SS=1;
for ii=1:length(F)
    
    
    realMemAccs(SS)=F(ii).realMemAcc;
    
    percSacColChgs(SS)=F(ii).percSacColChg;
    percSacOthers(SS)=F(ii).percSacOther;
    percSacNonObjs(SS)=F(ii).percSacNonObj;
    percSacUnspecs(SS)=F(ii).percSacUnspec;
    
    percFirstSacColChgs(SS)=F(ii).percFirstSacColChg;
    percFirstSacOthers(SS)=F(ii).percFirstSacOther;
    percFirstSacNonObjs(SS)=F(ii).percFirstSacNonObj;
    percFirstSacUnspecs(SS)=F(ii).percFirstSacUnspec;
    
    
    avgMemRTs(SS)=F(ii).avgMemRT;
    SS=SS+1;
    
end
end






