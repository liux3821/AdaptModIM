
%% HouseKeeping
clear all; close all;
commandwindow;
smallwin = input('Small Window y/n: ', 's');
testcond = input('Adaptor Same Orientation y/n: ', 's');
subject_code = input('Subject Code: ', 's');
exp_Loc = input('p for PsyPhys, l for Laptop: ', 's');
HideCursor;

%% Read code and save as CHAR F
%matlab saves F as numbers; needs char(F) to read;
%fid = fopen([mfilename '.m'],'r');
%completeScript = fread(fid)'; %#ok<*NASGU>
%eval(['code.ortmatching = char(completeScript);'])
%fclose(fid);

%% Adaptor Variables
% ORIENTATION
gratAngTest = [20, -20]; %EXPERIMENTAL CONDITION
gratAngCtrl = [20,20]; %CONTROL CONDITION
numPhases = 4;

%ANGULAR SIZE
%ATTN: ON LAPTOP SHOULD USE SMALLER ANGULAR SIZE:4.
%IN PSYCHOPHYSICS ROOM SHOULD USE ANGULAR SIZE 8.

adaptorsize = 4; %in degrees
testsize = 4; %in degrees

%SPATIAL FREQUENCY
%ON LAPTOP SHOULD USE CPD 4
%IN PSYCHOPHYSICS ROOM SHOULD USE CPD 1.5.

cpd = 1.5;
cpd_random = 3;
cps_adaptor = cpd*adaptorsize;
cps_test = cpd*testsize;
cps_random = cpd_random*adaptorsize;%Cycles per degree

%GABOR SIZE, MISCELLANEUOUS
gabWidth = 5;
testgabWidth = 150; %IMAGE GABOR SIZE
frameDur = 0.1;
testDur = 1; %IMAGE PRESENT TIME
maskSiz = 0.8;

%EXPERIMENT PARAMETERS
adapttime = 5; %specify preadaptation duration
topUpInTrial = 1;
topUpBTTrial = 3;
topup = 3;
fixpix = 8; %size of the fixation cross
%numoftrials = 10;

randomStartFrame = 2;
randomFinishFrame = (adapttime./frameDur)-10;
randomFrame = randi([randomStartFrame, randomFinishFrame]); %RANDOM FRAME TO PRESENT MULTIPLE SF

%VECTORS STORE RESULTS
taskTime = []; %GET TIME WHEN SUBJECT CONFIRM ORT MATCHING
timeSF = []; %vector stores time WHEN multiple SF appears
tasktimeSF = []; %vector stores time subjects RT to multiple SF
TAE = [];
%BASELINE RESULTS
TAE_base = [];
taskTime_base = [];
%% Task Variable
gratAngRef = 0;
task_size = adaptorsize; %in degrees
reference_size = task_size; %in degrees
cpd = 1.5;
cps_task = cpd*task_size;
cps_reference = cpd*reference_size;%Cycles per degree
gabWidth = 5;
fixpix = 8;
maskSiz = 0.8;
dur_Block = 360; %DURATION PER BLOCK
dur_Baseline = 5;

%% SOUND CUE VARIABLES
sndInfo.numSnds = 10;                                     % number of sound frequencies used to communicate the change in orientation
f0              = 400;                                    % Frequency (in Hz) of the orientation change cue
sf              = 44100;                                  % sampling frequency of the sound
sDur            = 0.01; %floor(1000*0.35*timeinfo.stimDur)/1000; % duration (in sec) of each sound sndInfo.blip
pDur            = 0.01; %floor(1000*0.35*timeinfo.stimDur)/1000; % duration (in sec) of silence between sndInfo.blips
sndInfo.totDur  = sDur+pDur;                              % total duration (in sec) of sndInfo.blip + silence
t               = 2*pi.*(1/sf:1/sf:sDur);

rampDur         = 0.0025;
ramp            = sin(linspace(0,pi/2,floor(rampDur*sf)));
rampTemplate    = [ramp ones(1,length(t)-2*length(ramp)) fliplr(ramp)];

for currSnd = 1:sndInfo.numSnds
    allFreqs(currSnd)       = 2^((currSnd-1)/12)*f0;
    sndInfo.blip(currSnd,:) = [rampTemplate.*sin(allFreqs(currSnd)*t) zeros(1,sf*pDur)]; % generate the sound itself and append it with desired duration of silence
end

% Perform basic initialization of the sound driver:
InitializePsychSound;

% Open 1 audio channel
sndInfo.pahandle = PsychPortAudio('Open', [], [], 1, sf, 1);
runMode = 1;
PsychPortAudio('RunMode', sndInfo.pahandle, runMode);

timeCheck = zeros(sndInfo.numSnds,2);
% Fill the audio playback buffer with the audio data in 'sndInfo.blip' and play each sound once:
for currSnd = 1:sndInfo.numSnds
    
    sndInfo.blipBuffer(currSnd) =  PsychPortAudio('CreateBuffer', sndInfo.pahandle, sndInfo.blip(currSnd,:));
    
end

PsychPortAudio('UseSchedule', sndInfo.pahandle, 1);
for currSnd = 1:sndInfo.numSnds
    loopStart = GetSecs;
    
    %ONE LEGO PIECE
    poppop(currSnd,sndInfo)
    timeCheck(currSnd,1) = (GetSecs - loopStart);
    WaitSecs(sndInfo.totDur);
    
end
% figure;subplot(1,2,1); hist(timeCheck(:,1)); subplot(1,2,2); hist(timeCheck(:,2))


%% Call Screen
AssertOpenGL;
Screen('Preference', 'VisualDebugLevel', 1); %get rid of the welcoming screen
screens=Screen('Screens');
screenNumber=max(screens);
screenRect=Screen('Rect', screenNumber);

%MONITOR PARAMETERS
if exp_Loc == 'p'
    monitor.viewDist = 49;
    monitor.size(1) = 40.5;
    monitor.size(2) = 29.5;
    monitor.center = screenRect(3:4)./2;
    offsetpix = 90;
    extraOffLeft = 16;
else if exp_Loc == 'l'
        monitor.viewDist = 58; %cm
        monitor.center = screenRect(3:4)./2;
        monitor.size(1) = 29; %Also cm
        monitor.size(2)= monitor.size(1).*(monitor.center(2)./monitor.center(1));
        offsetpix = 95;
        extraOffLeft = 15;
    end
end

%CONVERT VISUAL ANGLES INTO PIXELS
dpp = 180./(monitor.center(1).*2); %this is degree/pixel adjustment

%stimulus sizes
adaptorsizepixel = visAng2xyNew(adaptorsize,0,monitor);
%testsizepixel = visAng2xyNew(testsize,0,monitor);
tasksizepixel = visAng2xyNew(task_size,0,monitor);
refsizepixel = visAng2xyNew(reference_size,0,monitor);

%offsets sizes
offsetpixleft = (-1.*tasksizepixel)+extraOffLeft; %offset task grating to left
offsetpixup = adaptorsizepixel-offsetpix; %95 pixels if on laptap
offsetpixdown = -1.*offsetpixup;

white=WhiteIndex(screenNumber);
black=BlackIndex(screenNumber);
gray=round((white+black)/2);

if gray == white
    gray=white / 2;
end

inc=white-gray; %100% MICHELSON CONTRAST
taskCtrst = 0.5.*inc; %TEST CTRST:50%
adaptCtrst = 0.9.*inc;%ADAPTOR CTRST:90%

%% Calibration
if exp_Loc == 'p'
    addpath([pwd,'/loadCalibration']);
    load calib_macMini_NEC_CRT-0-22-Nov-2016.mat
else if exp_Loc == 'l'
        addpath([pwd,'/calibrationResults']);
        load calibration-0-27-Oct-2016.mat;
    end
end

dacsize = 8;  %By hand for now
monitor.gamInv = gamInv;
monitor.maxcol = 2.^dacsize-1;

bgRGB = [.5 .5 .5];
[bgLMS, ldRGB, actualdirLMS, ldContrast] = ...
    computeColorDirs(bgRGB, RGB2LMS, [1 1 1]);

[bgLMS, lmRGB, actualdirLMS, lmContrast] = ...
    computeColorDirs(bgRGB, RGB2LMS, [-1 1 0]);

newcmap = rgb2cmapramp(ldRGB,bgRGB,1.000,256,monitor.gamInv);
newclut = newcmap./monitor.maxcol;

%% Open Screen & Blend Function
% first two coordinates denote window position, second sets is window
if smallwin == 'y'
    [w,allcords] = Screen('OpenWindow',screenNumber, gray,[0,0,400,400]);
else
    [w,allcords] = Screen('OpenWindow',screenNumber, gray);
end
%alpha blending: add a top layer of gaussian envelope
Screen(w,'BlendFunction',GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
[xCenter,yCenter] = RectCenter(allcords); %get the center coordinates
%Fixation point parameter
fixRect1 = SetRect(0,0,fixpix,fixpix);
fixRect1 = CenterRect(fixRect1,allcords);
fixRect2 = SetRect(0,0,fixpix./2,fixpix./2);
fixRect2 = CenterRect(fixRect2,allcords);
imgRect_Adpt_Mid = CenterRectOnPoint([0,0,adaptorsizepixel,adaptorsizepixel],xCenter,yCenter);
%imgRect_Test = CenterRectOnPoint([0,0,testsizepixel,testsizepixel],xCenter,yCenter);
imgRect_Adpt(1,:) = OffsetRect(imgRect_Adpt_Mid,0,offsetpixup);
imgRect_Adpt(2,:) = OffsetRect(imgRect_Adpt_Mid,0,offsetpixdown);
%imgRect_Task = CenterRectOnPoint([0,0,refsizepixel,refsizepixel],xCenter,yCenter);
imgRect_Task = OffsetRect(imgRect_Adpt_Mid,offsetpixleft,0);
%imgRect_Test = OffsetRect(imgRect_Test,0,testsizepixel);

%DETERMINE GRATING ORIENTATIONS
if testcond == 'y'
    gratAng = gratAngCtrl;
else
    gratAng = gratAngTest;
end

for updown = 1:2
    for i=1:numPhases
        phase=(i/numPhases)*2*pi;
        adaptorGrating = drawGrating(gratAng(updown),adaptorsizepixel,...
            cps_adaptor,0,gabWidth,maskSiz,phase);
        randomGrating = drawGrating(gratAng(updown),adaptorsizepixel,...
            cps_random,0,gabWidth,maskSiz,phase);
        RandomGratingtex(updown,i)=Screen('MakeTexture', w, round(gray+...
            adaptCtrst*randomGrating));
        tex(updown,i)=Screen('MakeTexture', w, round(gray+...
            adaptCtrst*adaptorGrating));
    end
end

%READ IMAGES
d = dir('/Users/liux3821/Dropbox/MATLAB/Engel_Lab/Adaptation_TAE/testimages/image_vert*');
numofimages = length(d);
% % Use realtime priority for better timing precision:
priorityLevel=MaxPriority(w);
Priority(priorityLevel);

%MAKE REFERENCE GRATING TEXTURE
refGrating = drawGrating(gratAngRef,refsizepixel,cps_reference,0,gabWidth...
    ,maskSiz,0);
ref_tex = Screen('MakeTexture', w, round(gray+taskCtrst*refGrating));


%% PROMPT SCREEN
StartBaseline = 0;
while ~StartBaseline
    DrawFormattedText(w,['\n\n\n\n\n\n\nPlease Fixate at the Fixation Point During The Experiment\n'...
        'Press SPACE to start experiment'],'center','center',white);
    Screen('FillOval', w, gray, fixRect1);
    Screen('FillOval', w, black, fixRect2);
    Screen('Flip',w);
    [KeyDown,time,keyCode] = KbCheck;
    if KeyDown == 1
        if keyCode(KbName('SPACE'))
            StartBaseline = 1;
            break;
        end
    end
end
%% START EXPERIMENT
%tic;
baseLineStart = GetSecs; baseLineEnd = GetSecs; %TIME STAMP TO DETERMINE BASELINE DURATION

%% BASELINE ORT MATCHING MEASUREMENT
while baseLineEnd - baseLineStart < dur_Baseline
    %Get Mouse Movement to adjust Ort
    SetMouse(xCenter,yCenter,w);
    [xMouse,yMouse,buttons] = GetMouse(w);
    breakLoop = 0;
    
    %RANDOMIZE TASK LOCATION
    refRandomLoc = randi([1,2]);
    if refRandomLoc == 1
        imgRect_Reference = OffsetRect(imgRect_Adpt_Mid,0,offsetpixup);;
    else
        imgRect_Reference = OffsetRect(imgRect_Adpt_Mid,0,offsetpixdown);
    end
    
    %RANDOMIZE TASK ORIENTATION +/- 25DEG
    taskRandomOrt = randi([-25,25]);
    
    %TASK GRATING
    taskGrating = drawGrating(taskRandomOrt,tasksizepixel,cps_task,0,gabWidth...
        ,maskSiz,0);
    task_tex = Screen('MakeTexture', w, round(gray+taskCtrst*taskGrating));
    ort_change = 0; delta_ort = 0; old_ort = 0;
    %TRIAL LOOP
    while breakLoop == 0
        timeStamp = GetSecs;
        %TASK PRESENTATION: 250MS
        while GetSecs - timeStamp <= 0.25 && breakLoop == 0
            [ort_change,delta_ort,old_ort] = ORTchange(taskRandomOrt,dpp,xCenter,w,ort_change,delta_ort,old_ort);
            %txtOrt = num2str(ort_change);
            %txtOldOrt = num2str(old_ort);
            %txtOrtChange = num2str(ort_change);
            poppop(delta_ort,sndInfo);
            task_grating_change = drawGrating(ort_change,tasksizepixel,...
                cps_task,0,gabWidth,maskSiz,0);
            task_tex_change = Screen('MakeTexture', w, round(gray+taskCtrst*task_grating_change));
            Screen('DrawTexture', w, task_tex_change,[],imgRect_Task);
            Screen('DrawTexture', w, ref_tex,[],imgRect_Reference);
            %DrawFormattedText(w,['\n\n\n\n\n\n\n\n\n\n\n\nCurrent Ort is at\n',txtOrt,'degrees\n\n\n\n Current OLDOrt is at\n',txtOldOrt,'degrees\n\n\n\n Current DELTAOrt is at\n',txtOrtChange,'degrees\n\n\n\n'],'center','center',white);
            Screen('FillOval', w, gray, fixRect1);
            Screen('FillOval', w, black, fixRect2);
            Screen('Flip',w);
            %imageArray = Screen('GetImage',w);
            %imwrite(imageArray, 'ortmatching.jpg');
            [xMouse,yMouse,buttons] = GetMouse(w);
            [keyIsDown,foo,keyCode] = KbCheck;
            escapeKey(keyIsDown,keyCode);
            if buttons(1) == 1
                TAE_base(end+1) = [ort_change];
                taskTime_base(end+1) = GetSecs-baseLineStart;
                SetMouse(xCenter,yCenter,w);
                baseLineEnd = GetSecs;
                breakLoop = 1;
            end
        end
        
        timeStamp2 = GetSecs;
        %PRESENT BLANK SCREEN 200MS
        while GetSecs - timeStamp2 <= 0.2 && breakLoop == 0
            [ort_change,delta_ort,old_ort] = ORTchange(taskRandomOrt,dpp,xCenter,w,ort_change,delta_ort,old_ort);
            poppop(delta_ort,sndInfo);
            grayTex = Screen('MakeTexture',w,gray);
            Screen('DrawTexture',w,grayTex);
            Screen('FillOval', w, gray, fixRect1);
            Screen('FillOval', w, black, fixRect2);
            Screen('Flip',w);
            [xMouse,yMouse,buttons] = GetMouse(w);
            [keyIsDown,foo,keyCode] = KbCheck;
            escapeKey(keyIsDown,keyCode);
            if buttons(1) == 1
                TAE_base(end+1) = [ort_change];
                SetMouse(xCenter,yCenter,w);
                taskTime_base(end+1) = GetSecs-baseLineStart;
                baseLineEnd = GetSecs;
                breakLoop = 1;
            end
        end
    end
    Screen('DrawTexture',w,grayTex);
    Screen('FillOval', w, gray, fixRect1);
    Screen('FillOval', w, black, fixRect2);
    Screen('Flip',w);
    WaitSecs(2);
end

%% PROMPT SCREEN 2 TO PROCEED TO ADAPTATION
StartAdapt = 0;
while ~StartAdapt
    DrawFormattedText(w,['\n\n\n\n\n\n\nPlease Fixate at the Fixation Point During The Experiment\n'...
        'Press SPACE to continue'],'center','center',white);
    Screen('FillOval', w, gray, fixRect1);
    Screen('FillOval', w, black, fixRect2);
    Screen('Flip',w);
    [KeyDown,time,keyCode] = KbCheck;
    if KeyDown == 1
        if keyCode(KbName('SPACE'))
            StartAdapt = 1;
            break;
        end
    end
end

%% ADAPTATION
expStartTimeStamp = GetSecs; expEndTimeStamp = GetSecs; %TIME STAMPS TO DETERMINE BLOCK DURATION
trial_num = 1;

while expEndTimeStamp - expStartTimeStamp < dur_Block
    %Starting Parameters
    if trial_num == 1
        movieDurationSecs = adapttime;
    else
        movieDurationSecs = topUpBTTrial;
    end
    
    movieDurationFrames=round(movieDurationSecs./frameDur);
    movieFrameIndices = randi(numPhases,movieDurationFrames,1);
    
    for i=1:movieDurationFrames
        for updown = 1:2
            if i == randomFrame
                %Screen('DrawTexture', w, tex(updown, movieFrameIndices(i)),[],imgRect_Adpt(updown,:));
                Screen('DrawTexture',w,RandomGratingtex(updown, movieFrameIndices(i))...
                    ,[],imgRect_Adpt(updown,:));
            else
                Screen('DrawTexture', w, tex(updown, movieFrameIndices(i)),[],imgRect_Adpt(updown,:));
            end
        end
        Screen('FillOval', w, gray, fixRect1);
        Screen('FillOval', w, black, fixRect2);
        Screen('Flip', w);
        
        %GET MULTIPLE SF PRESENT TIME
        if i == randomFrame
            timeSF(end+1) = GetSecs;
        end
        WaitSecs(0.1);
        
        %CHECK KB DURING ADAPTATION
        [SFtask, time,keyCode] = KbCheck;
        if (SFtask == 1)
            if keyCode(KbName('SPACE'))
                tasktimeSF(end+1) = GetSecs; %GET TIME WHEN SUBJECT INDICATE MULTIPLE SF
            end
        end
        
        [keyIsDown,foo,keyCode] = KbCheck;
        escapeKey(keyIsDown,keyCode);
    end
    %% Start Task
    %Get Mouse Movement to adjust Ort
    SetMouse(xCenter,yCenter,w);
    [xMouse,yMouse,buttons] = GetMouse(w);
    breakLoop = 0;
    
    %RANDOMIZE TASK LOCATION
    refRandomLoc = randi([1,2]);
    if refRandomLoc == 1
        imgRect_Reference = OffsetRect(imgRect_Adpt_Mid,0,offsetpixup);;
    else
        imgRect_Reference = OffsetRect(imgRect_Adpt_Mid,0,offsetpixdown);
    end
    
    %RANDOMIZE TASK ORIENTATION +/- 25DEG
    taskRandomOrt = randi([-25,25]);
    
    %TASK GRATING
    taskGrating = drawGrating(taskRandomOrt,tasksizepixel,cps_task,0,gabWidth...
        ,maskSiz,0);
    task_tex = Screen('MakeTexture', w, round(gray+taskCtrst*taskGrating));
    
    
    %TRIAL LOOP
    while breakLoop == 0
        timeStamp = GetSecs;
        while GetSecs - timeStamp <= 0.25 && breakLoop == 0 && GetSecs - expEndTimeStamp > 0.25
            [ort_change,delta_ort,old_ort] = ORTchange(taskRandomOrt,dpp,xCenter,w,ort_change,delta_ort,old_ort);
            poppop(delta_ort,sndInfo);
            task_grating_change = drawGrating(ort_change,tasksizepixel,...
                cps_task,0,gabWidth,maskSiz,0);
            task_tex_change = Screen('MakeTexture', w, round(gray+taskCtrst*task_grating_change));
            Screen('DrawTexture', w, task_tex_change,[],imgRect_Task);
            Screen('DrawTexture', w, ref_tex,[],imgRect_Reference);
            Screen('FillOval', w, gray, fixRect1);
            Screen('FillOval', w, black, fixRect2);
            Screen('Flip',w);
            [keyIsDown,foo,keyCode] = KbCheck;
            escapeKey(keyIsDown,keyCode);
            [xMouse,yMouse,buttons] = GetMouse(w);
            if buttons(1) == 1
                TAE(end+1) = [ort_change];
                taskTime(end+1) = [GetSecs-expStartTimeStamp];
                SetMouse(xCenter,yCenter,w);
                expEndTimeStamp = GetSecs;
                trial_num = trial_num+1;
                breakLoop = 1;
            end
        end
        
        timeStamp2 = GetSecs;
        while GetSecs - timeStamp2 <= 0.2 && breakLoop == 0 &&...
                GetSecs - expEndTimeStamp > 0.25
            [ort_change,delta_ort,old_ort] = ORTchange(taskRandomOrt,dpp,xCenter,w,ort_change,delta_ort,old_ort);
            poppop(delta_ort,sndInfo);
            task_grating_change = drawGrating(ort_change,tasksizepixel,...
                cps_task,0,gabWidth,maskSiz,0);
            task_tex_change = Screen('MakeTexture', w, round(gray+taskCtrst*task_grating_change));
            grayTex = Screen('MakeTexture',w,gray);
            Screen('DrawTexture',w,grayTex);
            Screen('FillOval', w, gray, fixRect1);
            Screen('FillOval', w, black, fixRect2);
            Screen('Flip',w);
            [keyIsDown,foo,keyCode] = KbCheck;
            escapeKey(keyIsDown,keyCode);
            [xMouse,yMouse,buttons] = GetMouse(w);
            if buttons(1) == 1
                TAE(end+1) = [ort_change];
                taskTime(end+1) = [GetSecs-expStartTimeStamp];
                SetMouse(xCenter,yCenter,w);
                trial_num = trial_num+1;
                expEndTimeStamp = GetSecs;
                breakLoop = 1;
            end
        end
        
        %TOPUP: 2000MS
        if breakLoop == 0
            movieDurationFrames=round(topUpInTrial./frameDur);
            movieFrameIndices = randi(numPhases,movieDurationFrames,1);
            for i=1:movieDurationFrames
                for updown = 1:2
                    Screen('DrawTexture', w, tex(updown, movieFrameIndices(i)),[],imgRect_Adpt(updown,:));
                end
                Screen('FillOval', w, gray, fixRect1);
                Screen('FillOval', w, black, fixRect2);
                Screen('Flip', w);
                %WaitSecs(0.1);
                [ort_change,delta_ort,old_ort] = ORTchange(taskRandomOrt,dpp,xCenter,w,ort_change,delta_ort,old_ort);
                poppop(delta_ort,sndInfo);
                if buttons(1) == 1
                    TAE(end+1) = [ort_change];
                    taskTime(end+1) = [GetSecs-expStartTimeStamp];
                    SetMouse(xCenter,yCenter,w);
                    trial_num = trial_num+1;
                    expEndTimeStamp = GetSecs;
                    task_grating_change = drawGrating(ort_change,tasksizepixel,...
                        cps_task,0,gabWidth,maskSiz,0);
                    task_tex_change = Screen('MakeTexture', w, round(gray+taskCtrst*task_grating_change));
                    breakLoop = 1;
                end
                [keyIsDown,foo,keyCode] = KbCheck;
                escapeKey(keyIsDown,keyCode);
            end
        end
    end
    
    Screen('DrawTexture',w,grayTex);
    Screen('FillOval', w, gray, fixRect1);
    Screen('FillOval', w, black, fixRect2);
    Screen('Flip',w);
    WaitSecs(0.25);
end
%toc;

%% Save test result
resultFilePath = strcat(pwd,'/expResults/',subject_code,'.mat');
save(resultFilePath);

%% Close audio
PsychPortAudio('DeleteBuffer');
PsychPortAudio('Stop', sndInfo.pahandle);   % Stop audio playback
PsychPortAudio('Close', sndInfo.pahandle);  % Close the audio device

%% Close Screens
Screen('Close');
RestoreCluts;
Screen('CloseAll'); ShowCursor; ListenChar(0);
