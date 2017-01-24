function ascreaderGazeDirect(file,charsubj)
% function where inputs are filename of asc file and number of
% characters in the subject #, & to save workspace to .mat - JI
%modified from ascreaderNC_Inst
%% Find out which file and number of rows in file
nTrials=24; %fix this to 240 or what the new number is
center=[512 384];
memposcenter=[512 584;512 184; 338.8 484; 685.2 484; 685.2 284; 338.8 284];
nonObjCenter=[312 384; 712 384; 412 557.2; 412 210.8; 612 557.2; 612 210.8];
crit_ed=66;

r=0;
x=0;

subjnumdigit=charsubj;
subj=str2num(file(1:subjnumdigit));

% Open Data File

fid= fopen(file, 'rt');
% Loop through data file until we get a -1 indicating end of file
while(x~=(-1))
x=fgetl(fid);
r=r+1;
end
r = r-1;

disp(['Number of rows = ' num2str(r)])

%xy
objcenter=[1 2; 3 4; 5 6; 7 8; 9 10; 11 12];

%% PREANALYZE TO GET INFO ABOUT TRIAL
frewind(fid);
for i=1:r
name=fgetl(fid);
if regexpi(name, 'Trial=');
    trial=str2num(name(regexpi(name, 'Trial=\d', 'end'):size(name,2)));
end
if regexpi(name, 'evcodes');
    evco=regexpi(name, ' ', 'split');
    colorchange(trial)= str2num(evco{3})-100; %#ok<*SAGROW,*ST2NM>
    setsize(trial)= str2num(evco{4})-110;
    gazeDirect(trial)= str2num(evco{5})-120;
    whichSquare(trial)= str2num(evco{6})-130;
    memacc(trial)= str2num(evco{7})-140;
    memRT(trial)= str2num(evco{8});
    position1(trial)= str2num(evco{9})-150;
    position2(trial)= str2num(evco{10})-160;
    position3(trial)= str2num(evco{11})-170;
    position4(trial)= str2num(evco{12})-180;
    position5(trial)= str2num(evco{13})-190;
    position6(trial)=str2num(evco{14})-200;
    color1(trial)=str2num(evco{15})-210;
    color2(trial)=str2num(evco{16})-220;
    color3(trial)=str2num(evco{17})-230;
    color4(trial)=str2num(evco{18})-240;
    color5(trial)=str2num(evco{19})-250;
    color6(trial)=str2num(evco{20})-260;
    gazDirSeq1(trial)=str2num(evco{21})-270;
    gazDirSeq2(trial)=str2num(evco{22})-280;
    gazDirSeq3(trial)=str2num(evco{23})-290;
    gazDirSeq4(trial)=str2num(evco{24})-300;
    gazDirSeq5(trial)=str2num(evco{25})-310;
    gazDirSeq6(trial)=str2num(evco{26})-320;
    Obj_v_Non(trial)=str2num(evco{27})-330;
    gazePos1(trial)=str2num(evco{28})-340;
    gazePos2(trial)=str2num(evco{29})-350;
    gazePos3(trial)=str2num(evco{30})-360;
    gazePos4(trial)=str2num(evco{31})-370;
    gazePos5(trial)=str2num(evco{32})-380;
    gazePos6(trial)=str2num(evco{33})-390;
end
end


%% ANALYZE TRIALS FOR FIXATIONS
%Set up variables that I will write for each sample
fileName = ['ProcessedData/Gaze' num2str(subj)];
dataFile = fopen(fileName, 'a');
    header='Subject\tBlock\tcolChange\tTrial\tsetSiz\tgazDir\twhichSq\tmemAcc\tmemRT\tposition1\tp2\tp3\tp4\tp5\tp6\tcol\tcol2\tcol3\tcol4\tcol5\tcol6\tObjN\tGazPo1\tgp2\tgp3\tgp4\tgp5\tgp6\n';
    fprintf(dataFile, header);

ttstart=nan(nTrials,1);
GazeStart=nan(nTrials,1);
firstsactimestamp=nan(nTrials,1);
saconset=zeros(nTrials,1);
checknext=0;
gazePortion=0;
sac=0;

frewind(fid);
tt=0;
setok=0;

for i = 1:r

name = fgetl(fid); 

if setok==0 && (length(name)>4) && (strcmp(name(1:3),'MSG')) && (~isempty(regexpi(name, 'Trial=1')))
    setok=1;
elseif setok==1;

    %Look for Start of trial
    if (length(name)>4) && (strcmp(name(1:5),'START'));
        tt=tt+1;
        ttstart(tt)=i;
        gazePortion=0;
        firstsac=0; 
        checknext=0;
    end

    % Check if we're in the gaze direct portion of the trial, if so, set gazePortion=1
    if (length(name)>2) && (strcmp(name(1:3), 'MSG'));
         if isempty(regexpi(name, 'btGaz1'))==0 ; %Gaze Direction period onset 
             gazePortion=1;
             sac=1;
             firstsac=0;
             numfix=0;
             start=textscan(name, '%s %d %s');
             GazeStart(tt,1)=start{2};
    
             if(tt<21)
                 figure(tt);
                 
                 title(num2str(tt))
                 hold on; % hold figure
                 axis([0 center(1)*2 -center(2)*2  0]); %set figure axes
                 for beep=1:6
                     plot(memposcenter(beep,1), -memposcenter(beep,2),'c+')
                     plot(nonObjCenter(beep,1), -nonObjCenter(beep,2), 'cx')
                 end
             end
         elseif isempty(regexpi(name, 'MemTest'))==0; %Gaze Direct period offset 
             gazePortion=0;
             if ((numfix==0) && sac~=0);
                  fprintf(dataFile, ['%d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t' ...
                '%d\t %d\t %d\t %d\t %d\t %d\t' ...
                '%d\t %d\t %d\t %d\t %d\t %d\t' ...
                '%d\t %d\t %d\t %d\t %d\t %d\t %d\t' ...
                '%d\t %d\t %d\t %d\t %d\t %d\t'
                ],...
                subj, tt, colorchange(tt), max([saconset(tt,sac) 1]), -1, -1, -1, setsize(tt), gazeDirect(tt), whichSquare(tt),...
                memacc(tt), memRT(tt), position1(tt), position2(tt), position3(tt), position4(tt), position5(tt), position6(tt), color1(tt), color2(tt), color3(tt), ...
                color4(tt), color5(tt),color6(tt), gazDirSeq1(tt), gazDirSeq2(tt),gazDirSeq3(tt),gazDirSeq4(tt),gazDirSeq5(tt),gazDirSeq6(tt),Obj_v_Non(tt),gazePos1(tt),gazePos2(tt),gazePos3(tt),gazePos4(tt),gazePos5(tt),gazePos6(tt));
             end

         end
    end

    %Look for first saccade on the trial
    if (length(name)>4) && (strcmp(name(1:5), 'SSACC') && (gazePortion==1) && (firstsac==0)) %Note: This should be qualified with a test of a blink in the next samples
        sac=sac+1;
        sac1=textscan(name, '%s %s %d');
        saconset(tt,sac)= sac1{3}-GazeStart(tt);
        firstsactimestamp(tt,1)=sac1{3};
        firstsac=1;

    elseif (length(name)>4) && (strcmp(name(1:5), 'SSACC') && (gazePortion==1) && (firstsac==1))
        sac=sac+1;
        sac1=textscan(name, '%s %s %d');
        saconset(tt,sac)= sac1{3}-GazeStart(tt);
    end

    % if there is a blink start within 50ms of the first saccade onset, it is a
    % blink, and so re-set the calculation for first saccade onset
    if (length(name)>5) && (strcmp(name(1:6), 'SBLINK'))
        blink=textscan(name, '%s %s %d');
        if (blink{3}-firstsactimestamp(tt,1))< 50;
            firstsac=0;
            if sac>0
            sac=sac-1;
            end
        end
    end

    %Plot the location of the eyes on each trial, only plotting the sample
    if (length(name)>2) && (tt<21) && (~isempty(str2num(name(1:5)))) && (gazePortion==1) && (firstsac==1)
        
        currentsample=textscan(name, '%d %f %f %f %s');
        scatter(currentsample{2}, -currentsample{3});
    end


    %Look for fixations during gaze direct period, only after the first saccade occurred 
    if (length(name)>3) && (strcmp(name(1:4), 'SFIX') && (gazePortion==1) && (firstsac==1))
        numfix=numfix+1;
        checknext=i+1; %Now that fixation has started, pull out next sample to check location
    
    end
    
    if checknext==i;
        checknext=0;
        %pull out x and y values
        fix=textscan(name, '%d %d %d %d %s');
        eye(numfix,1,tt)=fix{2};
        eye(numfix,2,tt)=fix{3};
        %display([fix{2} ' ' fix{3}]);
        % if the fixation is still in the central fixation region, don't count
        % it and move on to the next fixation
        if (sqrt( double (( eye(numfix,1,tt) - center(1) )^2 + ( eye(numfix,2,tt) - center(2) )^2 )))<crit_ed;
            numfix=numfix-1;
        else
            %Compute the fixation location based on euclidean distance, using
            %specific critical distance, set to nan if not on an object
            fixloc(numfix,tt)=NaN;
            for loc=1:6
                edist(numfix,loc,tt)=(sqrt( double (( eye(numfix,1,tt) - memposcenter(loc,1) )^2 + ( eye(numfix,2,tt) - memposcenter(loc,2) )^2 )));
                edistNonObj(numfix,loc,tt)=(sqrt( double (( eye(numfix,1,tt) - nonObjCenter(loc,1) )^2 + ( eye(numfix,2,tt) - nonObjCenter(loc,2) )^2 )));

                if edist(numfix,loc,tt)<crit_ed
                    fixloc(numfix,tt)=loc;
                elseif edistNonObj(numfix,loc,tt)<crit_ed
                    fixloc(numfix,tt)=loc+6;
                end
            end
            
            %Determine if the fixation was on the color change object, other object, a
            %non-object gaze direct location, or no object
            %XX Critical Stuff!  
            if isnan(fixloc(numfix,tt))
                fixcode(numfix,tt)=4; % (other unspecified fixation location)
                nonclassifiedfix(numfix,tt)=min(edist(numfix,:,tt));
            elseif (isnan(fixloc(numfix,tt))==0) && (fixloc(numfix,tt)==whichSquare(tt))
                fixcode(numfix,tt)=1; % fixation at location of square that changes color
            elseif (isnan(fixloc(numfix,tt))==0) && (fixloc(numfix,tt)>6)
                fixcode(numfix,tt)=3; % fixation at non-object gaze direct location
            elseif isnan(fixloc(numfix,tt))==0
                fixcode(numfix,tt)=2; %fixation at any other of the 6 possible object locations
            end
            
            %header='Subject\tTrial\tColorChange\tSaccadeOnset\tFixNum\tFixCode\tFixLoc\tsetSize\tgazeDirect\twhichSquare\tmemAcc\tmemRT\tposition1\tp2\tp3\tp4\tp5\tp6\tcol\tcol2\tcol3\tcol4\tcol5\tcol6\tgazSeq1\tgSeq2\tgSeq3\tgSeq4\tgSeq5\tgSeq6\n';
            fprintf(dataFile, ['%d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t' ...
                '%d\t %d\t %d\t %d\t %d\t %d\t' ...
                '%d\t %d\t %d\t %d\t %d\t %d\t' ...
                '%d\t %d\t %d\t %d\t %d\t %d\t %d\t' ...
                '%d\t %d\t %d\t %d\t %d\t %d\t'
                ],...
                subj, tt, colorchange(tt), max([saconset(tt,sac) 1]), numfix, fixcode(numfix,tt), fixloc(numfix,tt), setsize(tt), gazeDirect(tt), whichSquare(tt),...
                memacc(tt), memRT(tt), position1(tt), position2(tt), position3(tt), position4(tt), position5(tt), position6(tt), color1(tt), color2(tt), color3(tt), ...
                color4(tt), color5(tt),color6(tt), gazDirSeq1(tt), gazDirSeq2(tt),gazDirSeq3(tt),gazDirSeq4(tt),gazDirSeq5(tt),gazDirSeq6(tt),Obj_v_Non(tt),gazePos1(tt),gazePos2(tt),gazePos3(tt),gazePos4(tt),gazePos5(tt),gazePos6(tt));
        end
    end
end

  
end
save(['ProcessedData/' file(1:2) '_workspace.mat'])

