%function GazeDirect
%GazeDirect (ColorVWMSac with Gaze Direction in ISI)
%--------------------------------------------------------------------------
% Written by Julianna Ianni Dec 2011
%  adapted from TVAEYE (N. Carlisle)
% uses function InitEyelinkDefaultsNC

try
    clear all;          % clear all matlab vars
    Screen('CloseAll'); % close any open psychtoolbox windows
    
    %%  PROMPT EXPERIMENTER
    
    promptTitle         = 'Experimental Setup Information';
    prompt              = {'Enter subject number: ', 'Block', 'Output'};
    promptNumAnsLines   = 1;
    promptDefaultAns    = {'99','1', 'RawData/'};
    answer              = inputdlg(prompt, promptTitle, promptNumAnsLines, promptDefaultAns);
    [expt.subjNum, expt.block, expt.saveDir] = deal(answer{:});
    
    subjNum = str2num(expt.subjNum);
    block = str2num(expt.block);
    ListenChar(2);
    % OPEN FILE AND PRINT TIME AND HEADER
    fileName = [expt.saveDir 'GazeDirect' expt.subjNum];
    dataFile = fopen(fileName, 'a');
    startTime = datestr(now, 0);
    fprintf(dataFile, '%s\t', startTime);
    fprintf(dataFile, '%i\n', subjNum);
    header='Subject\tBlock\tcolChange\tTrial\tsetSiz\tgazDir\twhichSq\tmemAcc\tmemRT\tpos1\tp2\tp3\tp4\tp5\tp6\tcol1\tc2\tc3\tc4\tc5\tc6\tSeq1\ts2\ts3\ts4\ts5\ts6\tObjN\tGazPo1\tgp2\tgp3\tgp4\tgp5\tgp6\n'
    fprintf(dataFile, header);
    
    %SETUP SCREEN AND TIMING CORRECTION
    bgColour=[145 145 145];
    oldRes=SetResolution(0,1024,768);
    [w, centerPt]=Screen('OpenWindow', 0, bgColour);
    expt.winPtr = w;
    centerx=centerPt(3)/2;
    centery=centerPt(4)/2;
    priorityLevel= MaxPriority(w);
    Priority(priorityLevel);
    timing=Screen('GetFlipInterval',w);% get the flip rate of current monitor.
    a='Resolution & Refresh ok';
    if (centerx~=512) || (timing>.012) || (timing<.011) %make sure we're running in 1024*768 at 85 Hz, else stop
        clear a;
        display('Please change screen to 1024x768 and 85Hz');
    end
    display a;
    
    timingcorrection=timing/2; %this ensures proper timing of the flips by making sure the command has enough time to execute before the next screen refresh
    [oldFontName,oldFontNumber]=Screen('TextFont', w, 'Helvetica');
    oldTextSize=Screen('TextSize', w, 12);
    
    %% SET UP EXPERIMENTAL VARIABLES
    %-------------------------------------------------------------------------
    numTrial=240% must be divisible by 2 and 3 currently
    colourIndex=[1 2 3 4 5 6 7];
    positionIndex=1:6;
    breaks=30; %take breaks after this many trials
    
    % set up object & array
    SquareSide=20;
    fixsize   = 10;
    fixwidth  = 4;
    fixspec= [centerx-.5*fixsize, centerx+.5*fixsize, centerx, centerx;   centery, centery, centery-.5*fixsize, centery+.5*fixsize];
    memradius=200;
    diag1=(memradius*(cos(pi/6)));
    diag2=(memradius*sin(pi/6));
    diagNonObj2=(memradius*cos(pi/6));
    diagNonObj1=(memradius*sin(pi/6));
    colour=[255 0 0; 0 255 0; 0 0 255; 255 0 255; 255 255 0; 0 0 0; 255 255 255]; %red, green, blue, violet, yellow, black, white
    fixcol=[0 0 0];
    possibleLocations=[centerx-SquareSide centery+memradius-SquareSide centerx+SquareSide centery+memradius+SquareSide; centerx-SquareSide centery-memradius-SquareSide centerx+SquareSide centery-memradius+SquareSide; centerx-diag1-SquareSide centery+diag2-SquareSide centerx-diag1+SquareSide centery+diag2+SquareSide; centerx+diag1-SquareSide centery+diag2-SquareSide centerx+diag1+SquareSide centery+diag2+SquareSide; centerx+diag1-SquareSide centery-diag2-SquareSide centerx+diag1+SquareSide centery-diag2+SquareSide; centerx-diag1-SquareSide centery-diag2-SquareSide centerx-diag1+SquareSide centery-diag2+SquareSide];
    nonObjLocs=[centerx-memradius-SquareSide centery-SquareSide centerx-memradius+SquareSide centery+SquareSide; centerx+memradius-SquareSide centery-SquareSide centerx+memradius+SquareSide centery+SquareSide; centerx-diagNonObj1-SquareSide centery+diagNonObj2-SquareSide centerx-diagNonObj1+SquareSide centery+diagNonObj2+SquareSide; centerx-diagNonObj1-SquareSide centery-diagNonObj2-SquareSide centerx-diagNonObj1+SquareSide centery-diagNonObj2+SquareSide; centerx+diagNonObj1-SquareSide centery+diagNonObj2-SquareSide centerx+diagNonObj1+SquareSide centery+diagNonObj2+SquareSide; centerx+diagNonObj1-SquareSide centery-diagNonObj2-SquareSide centerx+diagNonObj1+SquareSide centery-diagNonObj2+SquareSide];
    % set up times (note each is timing corrected to be called half a flip
    % interval before a flip)
    ITIbase     =1;
    ITIadd      =.4;
    ITI=zeros(numTrial);
    for a=1:numTrial;
        ITI(a)=round(((rand * ITIadd) + ITIbase)/timing)*timing-.5*timing;
    end
    fixdur =round(.5/timing)*timing-.5*timing;%time each fixation point should be presented for
    memdur =round(.5/timing)*timing-.5*timing;
    ISI=round(5/timing)*timing-.5*timing; %length of ISI
    maxtestdur = round(3/timing)*timing-.5*timing;   %time memory test is displayed
    artsupdur = round(1.5/timing)*timing-.5*timing;
    gazDirDur=round(.5/timing)*timing-.5*timing; %time each gaze direct cue is displayed
    btGazDur=round(((ISI-(6*gazDirDur))/7)/timing)*timing-.5*timing; %calculates "leftover" time in the ISI and makes evenly spaced gaps between the six gaze direct times
    
    
    % set up experimental variables
    %basic variables for non-catch trials
    colorchange=Shuffle([zeros(1,numTrial/2),ones(1,numTrial/2)]);
    setsize=Shuffle(repmat([1,3,6],1,numTrial/3));
    gazeDirect=Shuffle([zeros(1,numTrial/2),ones(1,numTrial/2)]);
    Obj_v_Non=Shuffle([zeros(1,numTrial/2),ones(1,numTrial/2)]); %whether gaze direct will be in object locations (1) or non-object locations (0);

    
    gazDirSeq=zeros(6,numTrial);
    colors=zeros(7,numTrial);
    positions=zeros(6,numTrial);
    gazePos=zeros(6,numTrial);
    colorset=zeros(3,6,numTrial);
    locationset=zeros(4,6,numTrial);
    gazeLocs=zeros(4,6,numTrial);
    gazeColorSet=zeros(3,6,numTrial);
    gazeColor=[30 30 30];
    whichSquare=zeros(1,numTrial);
    testcolorset=zeros(3,6,numTrial);
    for ii=1:numTrial
        colors(:,ii)=Shuffle(colourIndex);
        positions(:,ii)=Shuffle(positionIndex);
        gazDirSeq(:,ii)=Shuffle(1:6)';  %note that this does not refer to ABSOLUTE color or location; refers to positions and colors for a specific trial (used as index for gazeColorSet and locationset)
        if setsize(ii)==1
            colorset(:,:,ii)=[colour(colors(1,ii),:);bgColour;bgColour;bgColour;bgColour;bgColour]';
            if gazeDirect(ii)==1
                gazeColorSet(:,:,ii)=[gazeColor;bgColour;bgColour;bgColour;bgColour;bgColour]';
            else
                gazeColorSet(:,:,ii)=[bgColour;bgColour;bgColour;bgColour;bgColour;bgColour]';
            end
        elseif setsize(ii)==3
            colorset(:,:,ii)=[colour(colors(1,ii),:);colour(colors(2,ii),:);colour(colors(3,ii),:);bgColour;bgColour;bgColour]';
            if gazeDirect(ii)==1
                gazeColorSet(:,:,ii)=[gazeColor;gazeColor;gazeColor;bgColour;bgColour;bgColour]';
            else
                gazeColorSet(:,:,ii)=[bgColour;bgColour;bgColour;bgColour;bgColour;bgColour]';
            end
        else
            colorset(:,:,ii)=[colour(colors(1,ii),:);colour(colors(2,ii),:);colour(colors(3,ii),:);colour(colors(4,ii),:);colour(colors(5,ii),:);colour(colors(6,ii),:)]';
            if gazeDirect(ii)==1
                gazeColorSet(:,:,ii)=[gazeColor;gazeColor;gazeColor;gazeColor;gazeColor;gazeColor]';
            else
                gazeColorSet(:,:,ii)=[bgColour;bgColour;bgColour;bgColour;bgColour;bgColour]';
            end
        end
        locationset(:,:,ii)=[possibleLocations(positions(1,ii),:); possibleLocations(positions(2,ii),:); possibleLocations(positions(3,ii),:); possibleLocations(positions(4,ii),:); possibleLocations(positions(5,ii),:); possibleLocations(positions(6,ii),:)]';
        if Obj_v_Non(ii)==1 %set the gaze direct locations to object locations OR non object locations depending on trial type:
            gazeLocs(:,:,ii)=locationset(:,:,ii);
        else
            gazePos(:,ii)=Shuffle(positionIndex);
            gazeLocs(:,:,ii)=[nonObjLocs(gazePos(1,ii),:); nonObjLocs(gazePos(2,ii),:); nonObjLocs(gazePos(3,ii),:); nonObjLocs(gazePos(4,ii),:); nonObjLocs(gazePos(5,ii),:); nonObjLocs(gazePos(6,ii),:)]';
        end
        testcolorset(:,:,ii)=colorset(:,:,ii);
        if colorchange(ii)==1;
            if setsize(ii)==6
                choosing=Shuffle(1:6);
            elseif setsize(ii)==3
                choosing=Shuffle(1:3);
            else
                choosing=1;
            end
            whichSquare(ii)=choosing(1);
            testcolorset(:,whichSquare(ii),ii)=colour(colors(7,ii),:)';
        end
    end
    
    
    
    %% Output & Set Inital Variables
    save(['RawData/MFiles/Eye_' num2str(subjNum)]);
    
    memRT=zeros(1,numTrial);
    memRespCode=nan(1,numTrial);
    memacc=3*ones(1,numTrial);
    kcm=zeros(numTrial,256);
    kdm=zeros(1,numTrial);
    secm=zeros(1,numTrial);
    doublemessagetime=NaN(1,numTrial);
    
    %% Setup Eyelink System- must have eyelink computer started and running
    %%eyelink
    init=Eyelink('Initialize');
    el=EyelinkInitDefaultsNC(w);
    edffilename=[ expt.subjNum datestr(now,'mmddyy') '.edf'];
    openEDF=Eyelink('OpenFile', edffilename);
    
    HideCursor;
    
    result=EyelinkDoTrackerSetup(el,'c');
    DrawFormattedText(w, 'Which eye? 0=Left; 1=Right.', 'center', 'center');
    Screen('Flip',w);
    junk=NaN;
    while isnan(junk);
        [kde,sece,kce]=KbCheck;
        if (kde==1) && ((kce(48)==1) || (kce(49)==1) || (kce(96)==1) || (kce(97)==1))
            junk=1;
        elseif (kde==1) && (kce(81)==1)
            clear junk
        end
    end
    if ((kce(48)==1) || (kce(96)==1) ); eye=0;
    elseif ( (kce(49)==1) || (kce(97)==1));eye=1;
    else clear eye;
    end
    Screen('FillRect', w, bgColour);
    
    [vbl SOT]=Screen('Flip',w);
    if eye==0;
        Eyelink('Command', 'file_sample_filter=LEFT, GAZE, AREA, STATUS');
        Eyelink('Command', 'file_event_filter=LEFT,  FIXATION, SACCADE, BLINK, MESSAGE');
        Eyelink('Command', 'link_event_filter = LEFT, SACCADE,BLINK, MESSAGE');
    else
        Eyelink('Command', 'file_sample_filter=RIGHT, GAZE, AREA, STATUS');
        Eyelink('Command', 'file_event_filter= RIGHT, FIXATION, SACCADE, BLINK, MESSAGE');
        Eyelink('Command', 'link_event_filter = RIGHT,SACCADE,BLINK, MESSAGE');
    end
    
    screenpixstring=sprintf('screen_pixel_coords= %f,%f,%f,%f',centerPt);      %you don't want eyelink to take those from physical.ini, because that means changing physical.ini for every screen resolution setting
    Eyelink('Command', 'clear_screen 0');
    
    
    %% TRIALS
    %--------------------------------------------------------------------------
    for tt=1:numTrial  %Testing
        display(tt);
        
        %wait before starting the experiment-----
        if tt==1
            DrawFormattedText(w, 'Press the space to begin.', 'center', 'center');
            Screen('Flip',w);
            junk=NaN;
            beforestart=GetSecs;
            while isnan(junk) && GetSecs-beforestart<30;
                [kd,sec,kc]=KbCheck;
                if (kd==1) && (kc(32)==1)
                    junk=1;
                elseif (kd==1) && (kc(81)==1)
                    clear junk
                end
            end
            Screen('FillRect', w, bgColour);
            [vbl SOT]=Screen('Flip',w);
        end
        
        %Take a break every
        %"breaks" trials
        if (mod(tt,breaks)==1) && (tt~=1)
            
            DrawFormattedText(w, 'Press the space to begin the next set.', 'center', centery+150);
            Screen('Flip',w);
            junky=NaN;
            WaitSecs(2);
            takeabreak=GetSecs;
            while isnan(junky) && GetSecs-takeabreak<60;
                [kd,sec,kc]=KbCheck;
                if (kd==1) && (kc(32)==1)
                    junky=1;
                elseif (kd==1) && (kc(81)==1)
                    clear junky
                end
            end
            Screen('FillRect', w, bgColour);
            [vbl SOT]=Screen('Flip',w);
        end
        %Drift Correction
        drift=EyelinkDoDriftCorrect(el);
        
        %Start Recording Eye for current trial
        status=['record_status_message "Trial ' num2str(tt) '"'];
        Eyelink('Command', status);
        Eyelink('Message', '%s%d', 'Trial=', tt);
        
        Screen('FillRect', w, bgColour);
        [vbl SOT]=Screen('Flip',w);
        
        Eyelink('StartRecording');
        
        %Display Articulatory Suppression Stimuli
        DrawFormattedText(w, num2str(randsample(1:9,2)), 'center', 'center');
        [vbl SOT] = Screen('Flip',w,SOT+ITI(tt));
        
        %Display fixation-------------------
        Screen('FillRect', w, bgColour);
        Screen('DrawLines', w, fixspec, fixwidth, fixcol);
        [vbl SOT] = Screen('Flip',w,SOT+artsupdur);
        Eyelink('Message', 'initfix'); %InitialFixationCross
        
        %Display memory -----------
        Screen('FillRect',w, bgColour);
        Screen('DrawLines', w, fixspec, fixwidth, fixcol);
        
        Screen('FillRect',w,colorset(:,:,tt),locationset(:,:,tt));
        
        [vbl SOT]=Screen('Flip',w, SOT+fixdur);
        Eyelink('Message', 'memory'); %Memory array
        
        %Display fixation--------------------
        Screen('FillRect', w, bgColour);
        Screen('DrawLines', w, fixspec, fixwidth, fixcol);
        [vbl SOT] = Screen('Flip',w,SOT+memdur);
        Eyelink('Message','btGaz1'); %fixation before gaze direct 1
        
        
        %display gaze direct cues in a random sequence:
        
        %Display Gaze Direct 1 or blank---------------
        Screen('FillRect',w, bgColour);
        Screen('DrawLines', w, fixspec, fixwidth, fixcol);
       Screen('FrameRect', w, gazeColorSet(:,gazDirSeq(1,tt),tt)', (gazeLocs(:,gazDirSeq(1,tt),tt)'));
       [vbl SOT]= Screen('Flip',w,SOT+btGazDur);
        Eyelink('Message', 'gaz1'); %Gaze Direct 1
        
        %Display fixation--------------------
        Screen('FillRect', w, bgColour);
        Screen('DrawLines', w, fixspec, fixwidth, fixcol);
        %[vbl SOT] = Screen('Flip',w,SOT+gazDirDur);
        [vbl SOT] = Screen('Flip',w,SOT+gazDirDur);
        Eyelink('Message', 'btGaz2'); %fixation before gaze direct 2
        
        %Display Gaze Direct 2 or blank---------------
        Screen('FillRect',w, bgColour);
        Screen('DrawLines', w, fixspec, fixwidth, fixcol);
        Screen('FrameRect',w,gazeColorSet(:,gazDirSeq(2,tt),tt),gazeLocs(:,gazDirSeq(2,tt),tt));
        [vbl SOT]= Screen('Flip',w,SOT+btGazDur);
        Eyelink('Message', 'gaz2'); %Gaze Direct 2
        
        %Display fixation--------------------
        Screen('FillRect', w, bgColour);
        Screen('DrawLines', w, fixspec, fixwidth, fixcol);
        [vbl SOT] = Screen('Flip',w,SOT+gazDirDur);
        Eyelink('Message', 'btGaz3'); %fixation before gaze direct 3
        
        %Display Gaze Direct 3 or blank---------------
        Screen('FillRect',w, bgColour);
        Screen('DrawLines', w, fixspec, fixwidth, fixcol);
        Screen('FrameRect',w,gazeColorSet(:,gazDirSeq(3,tt),tt),gazeLocs(:,gazDirSeq(3,tt),tt));
        [vbl SOT]= Screen('Flip',w,SOT+btGazDur);
        Eyelink('Message', 'gaz3'); %Gaze Direct 3
        
        %Display fixation--------------------
        Screen('FillRect', w, bgColour);
        Screen('DrawLines', w, fixspec, fixwidth, fixcol);
        [vbl SOT] = Screen('Flip',w,SOT+gazDirDur);
        Eyelink('Message', 'btGaz4'); %fixation before gaze direct 4
        
        %Display Gaze Direct 4 or blank---------------
        Screen('FillRect',w, bgColour);
        Screen('DrawLines', w, fixspec, fixwidth, fixcol);
        Screen('FrameRect',w,gazeColorSet(:,gazDirSeq(4,tt),tt),gazeLocs(:,gazDirSeq(4,tt),tt));
        [vbl SOT]= Screen('Flip',w,SOT+btGazDur);
        Eyelink('Message', 'gaz4'); %Gaze Direct 4
        
        %Display fixation--------------------
        Screen('FillRect', w, bgColour);
        Screen('DrawLines', w, fixspec, fixwidth, fixcol);
        [vbl SOT] = Screen('Flip',w,SOT+gazDirDur);
        Eyelink('Message', 'btGaz5'); %fixation before gaze direct 5
        
        %Display Gaze Direct 5 or blank---------------
        Screen('FillRect',w, bgColour);
        Screen('DrawLines', w, fixspec, fixwidth, fixcol);
        Screen('FrameRect',w,gazeColorSet(:,gazDirSeq(5,tt),tt),gazeLocs(:,gazDirSeq(5,tt),tt));
        [vbl SOT]= Screen('Flip',w,SOT+btGazDur);
        Eyelink('Message', 'gaz5'); %Gaze Direct 5
        
        %Display fixation--------------------
        Screen('FillRect', w, bgColour);
        Screen('DrawLines', w, fixspec, fixwidth, fixcol);
        [vbl SOT] = Screen('Flip',w,SOT+gazDirDur);
        Eyelink('Message', 'btGaz6'); %fixation before gaze direct 6
        
        %Display Gaze Direct 6 or blank---------------
        Screen('FillRect',w, bgColour);
        Screen('DrawLines', w, fixspec, fixwidth, fixcol);
        Screen('FrameRect',w,gazeColorSet(:,gazDirSeq(6,tt),tt),gazeLocs(:,gazDirSeq(6,tt),tt));
        [vbl SOT]= Screen('Flip',w,SOT+btGazDur);
        Eyelink('Message', 'gaz6'); %Gaze Direct 6
        
        %Display fixation--------------------
        Screen('FillRect', w, bgColour);
        Screen('DrawLines', w, fixspec, fixwidth, fixcol);
        [vbl SOT] = Screen('Flip',w,SOT+gazDirDur);
        Eyelink('Message', 'btGaz7'); %fixation AFTER gaze direct 6
        
        %Display TEST ---------------------
        Screen('FillRect', w, bgColour);
        
        Screen('FillRect',w,testcolorset(:,:,tt),locationset(:,:,tt));
        Screen('DrawLines', w,fixspec , fixwidth, fixcol);
        [vbl SOT]= Screen('Flip',w,SOT+btGazDur);
        Eyelink('Message', 'MemTest'); %test
        
        %Record memory response-------------
        junk=nan;
        beforememresp=GetSecs;
        while isnan(junk) && GetSecs-beforememresp<maxtestdur;
            [kdm,secm,kcm]=KbCheck;
            if ((kcm(38)==1) || (kcm(40)==1)) %same/up or different/down
                junk=1;
            elseif (kd==1) && (kc(81)==1)
                clear junk
            end
        end
        Eyelink('Message', 'mresponse'); %Memory Response
        memRT(tt) = 1000*(secm-beforememresp);
        if kcm(38)== 1; memRespCode(tt)=0; %respond same
        elseif kcm(40)==1; memRespCode(tt)=1; %respond different
        else memRespCode(tt)=NaN;
        end
        if colorchange(tt)==memRespCode(tt);  memacc(tt)=1;
        elseif isnan(memRespCode(tt)); memacc(tt)=2; %if there's no response
        else memacc(tt)=0;
        end
        
        %Blank Screen ITI
        Screen('FillRect', w, bgColour);
        [vbl SOT]= Screen('Flip',w);
        
        
        
        
        
        doublemessagestart=GetSecs;
        Eyelink('Message', '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 'colChg', 'setSiz', 'gazDir', 'whiSq','Acc', 'RT', 'pos1', 'p2','p3', 'p4','p5','p6','c1', 'c2', 'c3','c4','c5','c6','Seq1','S2','S3','S4','S5','S6','ObjN','gazPo1','g2','g3','g4','g5','g6');
        Eyelink('Message', '%s %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d', 'Devcodes', 100+colorchange(tt), 110+setsize(tt), 120+gazeDirect(tt), 130+whichSquare(tt), 140+memacc(tt), round(memRT(tt)), 150+positions(1,tt), 160+positions(2,tt), 170+positions(3,tt), 180+positions(4,tt), 190+positions(5,tt), 200+positions(6,tt), 210+colors(1,tt), 220+colors(2,tt), 230+colors(3,tt), 240+colors(4,tt), 250+colors(5,tt), 260+colors(6,tt), 270+gazDirSeq(1,tt), 280+gazDirSeq(2,tt), 290+gazDirSeq(3,tt), 300+gazDirSeq(4,tt), 310+gazDirSeq(5,tt), 320+gazDirSeq(6,tt), 330+Obj_v_Non(tt),340+gazePos(1,tt),350+gazePos(2,tt),360+gazePos(3,tt),370+gazePos(4,tt),380+gazePos(5,tt),390+gazePos(6,tt));
        doublemessagetime(tt)=GetSecs-doublemessagestart;
        Eyelink('StopRecording');
        %header='Subject\tBlock\tcolChange\tTrial\tsetSiz\tgazDir\twhichSq\tmemAcc\tmemRT\tpos1\tp2\tp3\tp4\tp5\tp6\tcol1\tc2\tc3\tc4\tc5\tc6\tSeq1\ts2\ts3\ts4\ts5\ts6\tObjN\tGazPo1\tgp2\tgp3\tgp4\tgp5\tgp6\n'
        
        
        fprintf(dataFile, ['%d\t %d\t %d\t %d\t %d\t %d\t %d\t' ...
            '%d\t %d\t %d\t %d\t %d\t %d\t %d\t' ...
            '%d\t %d\t %d\t %d\t %d\t %d\t' ...
            '%d\t %d\t %d\t %d\t %d\t %d\t %d\t' ...
            '%d\t %d\t %d\t %d\t %d\t %d\t %d\n'], ...
            subjNum, block, colorchange(tt), tt, setsize(tt), gazeDirect(tt), whichSquare(tt), ...
            memacc(tt), memRT(tt),positions(1,tt),positions(2,tt),positions(3,tt),positions(4,tt),positions(5,tt),positions(6,tt), ...
            colors(1,tt),colors(2,tt),colors(3,tt),colors(4,tt),colors(5,tt),colors(6,tt), ...
            gazDirSeq(1,tt),gazDirSeq(2,tt),gazDirSeq(3,tt),gazDirSeq(4,tt),gazDirSeq(5,tt),gazDirSeq(6,tt), ...
            Obj_v_Non(tt),gazePos(1,tt),gazePos(2,tt),gazePos(3,tt),gazePos(4,tt),gazePos(5,tt),gazePos(6,tt));
        
        
    end
    DrawFormattedText(w, 'You''re done! Please call for the experimenter.','center', 'center');
    Screen('Flip',w);
    junk=NaN;
    WaitSecs(5);
    
    
    
    ShowCursor;
    ListenChar(1);
    priority(0);
    Screen('CloseAll');
    Eyelink('CloseFile');
    statusrec=Eyelink('ReceiveFile', edffilename);
    Eyelink('Shutdown')
    exptRes=SetResolution(0,1280,1024);
    endTime = datestr(now, 0);
    display([startTime,'  ', endTime]);
    
catch ME
    ShowCursor;
    ListenChar(1);
    priority(0);
    Screen('CloseAll')
    ME.stack
    ME.message
    Eyelink('CloseFile');
    statusrec=Eyelink('ReceiveFile', edffilename);
    Eyelink('Shutdown')
    exptRes=SetResolution(0,1280,1024);
    endTime = datestr(now, 0);
    display([startTime,'  ', endTime]);
end



