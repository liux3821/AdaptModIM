% To do:
% remember to load calibration mat

% Problems:
%{
% control adaptors with a different orientation -- pending
% contrast judging strategies: blackest || whitest || mean -- raise again;
control by gray scale
% steeper & narrower gaussian markers (Mark) -- solved by using small
gaussian and big mask
% contrast weber's law; not change much for high var (Juraj) -- solved by
using smaller bins
%}

%% Inputs
clear all; close all; commandwindow;
sbjName = input('Name: ', 's');
expCondArray = input(['Enter exp condition sequence (no seperation)',...
    '\n(1: high variance; 2: low variance): '],'s');
expStageArray = '132323';% '132323';

ListenChar(2); HideCursor;

%% Calibration
addpath([pwd,'/calibrationResults']);
load calibration-0-22-Nov-2016.mat;
dacsize = 8;
monitor.gamInv = gamInv;
monitor.maxcol = 2.^dacsize-1;

bgRGB = [.5 .5 .5];
[bgLMS, ldRGB, actualdirLMS, ldContrast] = ...
    computeColorDirs(bgRGB, RGB2LMS, [1 1 1]);

[bgLMS, lmRGB, actualdirLMS, lmContrast] = ...
    computeColorDirs(bgRGB, RGB2LMS, [-1 1 0]);

newcmap = rgb2cmapramp(ldRGB,bgRGB,1.000,256,monitor.gamInv);
newclut = newcmap./monitor.maxcol;

%% open window
% Screen('Preference', 'SkipSyncTests', 1);
grey = [128,128,128]; white = [255,255,255]; black = [0,0,0];
screenNumbers = Screen('Screens');
screenNumber = max(screenNumbers);% the external screen if not mirrored
[winPtr,winRect] = Screen('OpenWindow',screenNumber,grey);
[winCenterX,winCenterY] = RectCenter(winRect);

% Load color lookup table (clut)
oldclut = LoadIdentityClut(winPtr);
Screen('LoadNormalizedGammaTable',winPtr,newclut,2);

%% image size & location
% for CRT monitor
visDiscm = 49;
monWidcm = 29.5;
monWidpix = winRect(4);
monLencm = 40.5 ;
monLenpix = winRect(3);

imgVaAdpt = 5.4 / 180 * pi; % larger gratings for adptators
imgWidPixAdpt = round(tan(imgVaAdpt) * visDiscm / monWidcm * monWidpix);
imgLenPixAdpt = round(tan(imgVaAdpt) * visDiscm / monLencm * monLenpix);
% from Min (2012)
imgVaTest = 4.5 / 180 * pi;
imgWidPixTest = round(tan(imgVaTest) * visDiscm / monWidcm * monWidpix);
imgLenPixTest = round(tan(imgVaTest) * visDiscm / monLencm * monLenpix);
imgLocDevVa = 4.1 / 180 * pi;
imgLocDevPix = round(tan(imgLocDevVa) * visDiscm / monLencm * monLenpix);

testSide = 'left'; matchSide = 'right';
if testSide == 'left'
    imgRect_Adpt = CenterRectOnPoint([0,0,imgWidPixAdpt,imgLenPixAdpt],winCenterX-imgLocDevPix,winCenterY);
    imgRect_Test = CenterRectOnPoint([0,0,imgWidPixTest,imgLenPixTest],winCenterX-imgLocDevPix,winCenterY);
    imgRect_Match = CenterRectOnPoint([0,0,imgWidPixTest,imgLenPixTest],winCenterX+imgLocDevPix,winCenterY);
else
    imgRect_Adpt = CenterRectOnPoint([0,0,imgWidPixAdpt,imgLenPixAdpt],winCenterX+imgLocDevPix,winCenterY);
    imgRect_Test = CenterRectOnPoint([0,0,imgWidPixTest,imgLenPixTest],winCenterX+imgLocDevPix,winCenterY);
    imgRect_Match = CenterRectOnPoint([0,0,imgWidPixTest,imgLenPixTest],winCenterX-imgLocDevPix,winCenterY);
end

fixRect_w = CenterRectOnPoint([0,0,10,10],winCenterX,winCenterY);
fixRect_b = CenterRectOnPoint([0,0,4,4],winCenterX,winCenterY);


%% dur & control
nTrial_b = 180; % 3 min
nTrial_adpt = 60; % 60*4.9~=5min
nTrial_dadpt = 60; % 60*4.9~=5min

Dur_Adpt = 4; Dur_DynmAdpt = 0.2; nDynmAdpt = Dur_Adpt/Dur_DynmAdpt;
Dur_Blank1 = 0.35;
Dur_Test = 0.2;
Dur_Blank2 = 0.35;

flipAheadTime = 0.008;

%% distribution statistics
% for adaptors
unfpara_HV_A = [0.05,0.15,0.05;0.40,0.50,0.15;0.75,0.85,0.8];% barWid:0.1
unfpara_HV_D = [0.05,0.15,0.8;0.40,0.50,0.15;0.75,0.85,0.05];% barWid:0.1
unfpara_LV_A = [.7125-0.025,.7125+0.025,1];% mean: .7125, barWid: .05
unfpara_LV_D = [.1875-0.025,.1875+0.025,1];% mean: .1875, barWid: .05

% for test
ctrArray_Test = 0.25; 

% for match scale
ctrArray_Match=linspace(0,0.5,51); 

%% other grating properties
ort = pi/2;
sf = 1.5; % degree, from Min (2012)
cyclAdpt = sf * (imgVaAdpt*180/pi); % sf * va(in degree)
cyclTest = sf * (imgVaTest*180/pi); % sf * va(in degree)
bgc = 128;
gsStr = 0.1;
maskPrp = 0.8; % small gaussian and larger mask to mask only at edges

%% define keyboard
keyTest = KbName('LeftArrow'); keyPrsCodeTest = 1;
keyMatch = KbName('RightArrow'); keyPrsCodeMatch = 2;
keyEsc = KbName('ESCAPE'); keySpace = KbName('SPACE');

%% PREPARE Stimuli Matrix
rng('shuffle'); randseed = rng;

imgMtrAdpt = cell(1,nDynmAdpt);
for i = 1:nDynmAdpt
    imgMtrAdpt{i} = makeGrating(imgWidPixAdpt,ort,cyclAdpt,rand*2*pi,imgWidPixAdpt*maskPrp,gsStr,0);
    imgMtrTest{i} = makeGrating(imgWidPixTest,ort,cyclTest,rand*2*pi,imgWidPixTest*maskPrp,gsStr,0);
end
save imgMtr imgMtrAdpt imgMtrTest;

load imgMtr imgMtrAdpt; load imgMtr imgMtrTest; 
%% PRESENT Adaptors->Blank1->Test&Match->Blank2

for q = 1:length(expCondArray)
    % define rec and first ctrToPresent for each cond
    Rec = [];
    
    expCond = expCondArray(q);
    DrawFormattedText(winPtr,['Please fixate on the center point and'],'center',winCenterY-30,white);
    DrawFormattedText(winPtr,['Press ',KbName(keyTest),' if the image on the ',testSide,' has higher contrast.'],'center',winCenterY+0,white);
    DrawFormattedText(winPtr,['Press ',KbName(keyMatch),' if the image on the ',matchSide,' has higher contrast.'],'center',winCenterY+30,white);
    DrawFormattedText(winPtr,['Press SPACE to start.'],'center',winCenterY+60,white);
    flipStartTime = Screen('Flip',winPtr);
    while 1
        [keyIsDown,~,keyCode] = KbCheck(-3);
        if keyIsDown && keyCode(keySpace)
            break
        end
    end
    
    % a 1s blank after space key
    tempFlipStartTime = Screen('Flip',winPtr);
    Screen('Flip',winPtr,tempFlipStartTime+1);
    
    stiStart = Screen('Flip',winPtr);    
    for r = 1:length(expStageArray)                
        %% define contrast array for specific cond/stage
        % initial match contrast
        ctrMatchToPresent = ctrArray_Test;
        
        % adpt contrast
        expStage = expStageArray(r);
        if expCond == '1'
            if expStage == '1'
                ctrArray_Adpt = ones(nTrial_b,nDynmAdpt)*0;
            elseif expStage == '2'                
                ctrArray_Adpt = reshape(splFrmMultipleUnf(unfpara_HV_A,nTrial_adpt*nDynmAdpt),nTrial_adpt,nDynmAdpt);                
            elseif expStage == '3'                
                ctrArray_Adpt = reshape(splFrmMultipleUnf(unfpara_HV_D,nTrial_adpt*nDynmAdpt),nTrial_adpt,nDynmAdpt);
            end
        elseif expCond == '2'
            if expStage == '1'
                ctrArray_Adpt = ones(nTrial_b,nDynmAdpt)*0;
            elseif expStage == '2'                
                ctrArray_Adpt = reshape(splFrmMultipleUnf(unfpara_LV_A,nTrial_adpt*nDynmAdpt),nTrial_adpt,nDynmAdpt);                
            elseif expStage == '3'                
                ctrArray_Adpt = reshape(splFrmMultipleUnf(unfpara_LV_D,nTrial_adpt*nDynmAdpt),nTrial_adpt,nDynmAdpt);
            end      
        end
                
        % prepare keyboard record
        Rec_temp = [];
        isKeyRecordTest = 1;
        isKeyRecordMatch = 1;
        
        switch expStage
            case '1'
                testStartCode = 100;
            case '2'
                testStartCode = 200;
            case '3'
                testStartCode = 300;
        end                       
        
        for i = 1:size(ctrArray_Adpt,1)
            %% present Adaptors
            if expStage == '1'
                timeCtrl = 1/nDynmAdpt; % gap for baseline: 1s
            else
                timeCtrl = Dur_DynmAdpt; 
            end
            
            for j = 1:nDynmAdpt
                imgMtrUpdated = bgc + ((bgc-1)*ctrArray_Adpt(i,j)).*imgMtrAdpt{randi(nDynmAdpt,1)}; % map (-1,1) to a range centered on 128 (depending on ctr)
                imgTex_Adpt = Screen('MakeTexture',winPtr,imgMtrUpdated);
                Screen('DrawTexture',winPtr,imgTex_Adpt,[],imgRect_Adpt);
                Screen('FillOval',winPtr,white,fixRect_w);
                Screen('FillOval',winPtr,black,fixRect_b);
                
                % control Dur
                flipStartTime = Screen('Flip',winPtr);
                while 1
                    timeTmp = GetSecs;
                    if (timeTmp - flipStartTime) > timeCtrl - flipAheadTime
                        break;
                    end
                    [keyIsDown,keySecs,keyCode] = KbCheck(-3);
                    
                    if keyIsDown == 0
                        isKeyRecordTest = 1;
                        isKeyRecordMatch = 1;
                    end
                                        
                    if keyIsDown && keyCode(keyTest) && isKeyRecordTest
                        Rec_temp = [Rec_temp;keyPrsCodeTest,0,keySecs];
                        isKeyRecordTest = 0;
                    end
                    if keyIsDown && keyCode(keyMatch) && isKeyRecordMatch
                        Rec_temp = [Rec_temp;keyPrsCodeMatch,0,keySecs];
                        isKeyRecordMatch = 0;
                    end
                    if keyIsDown && keyCode(keyEsc)
                        Rec = [Rec; Rec_temp];
                        Rec(:,end) = Rec(:,end) - stiStart;                                    
                        temp = clock; save(['ctrKeyRec_Cond',expCond,'_',sbjName,'_',num2str(temp(1:end-1))],'Rec','randseed');
                        RestoreCluts; Screen('CloseAll'); ShowCursor;ListenChar(0);
                        return;
                    end
                end
            end
            %% present Blank1
            Screen('FillOval',winPtr,white,fixRect_w);
            Screen('FillOval',winPtr,black,fixRect_b);
            
            % control duration
            flipStartTime = Screen('Flip',winPtr);
            while 1
                timeTmp = GetSecs;
                if (timeTmp - flipStartTime) > Dur_Blank1 - flipAheadTime
                    break;
                end                
                [keyIsDown,keySecs,keyCode] = KbCheck(-3);
                
                if keyIsDown == 0
                    isKeyRecordTest = 1;
                    isKeyRecordMatch = 1;
                end
                                
                if keyIsDown && keyCode(keyTest) && isKeyRecordTest
                    Rec_temp = [Rec_temp;keyPrsCodeTest,0,keySecs];
                    isKeyRecordTest = 0;
                end
                if keyIsDown && keyCode(keyMatch) && isKeyRecordMatch
                    Rec_temp = [Rec_temp;keyPrsCodeMatch,0,keySecs];
                    isKeyRecordMatch = 0;
                end
                if keyIsDown && keyCode(keyEsc)
                    Rec = [Rec; Rec_temp];
                    Rec(:,end) = Rec(:,end) - stiStart;                                    
                    temp = clock; save(['ctrKeyRec_Cond',expCond,'_',sbjName,'_',num2str(temp(1:end-1))],'Rec','randseed');                    
                    RestoreCluts; Screen('CloseAll'); ShowCursor;ListenChar(0);
                    return;
                end
            end
            
            %% present Test&Match
            if isempty(Rec_temp) % if no keyRec
                J = 1; % assume choose target as of high ctr for the first trial
            else
                J = Rec_temp(max(find(ismember(Rec_temp(:,1),[keyPrsCodeTest,keyPrsCodeMatch]))),1); % find the latest keyboard input
                if isempty(J) % if no keyPress                    
                    J = randi(2); % randomly choose a response
                end
            end
            if i ~= 1
                nextCtrToPresent = rvs3StairCase(ctrArray_Test,ctrMatchToPresent,J,ctrArray_Match,5,5,3,3,1,1); % initial step: 5*0.02 = 0.1
                ctrMatchToPresent = [ctrMatchToPresent,nextCtrToPresent];
            else
                nextCtrToPresent = ctrArray_Test;
                ctrMatchToPresent = [ctrMatchToPresent,nextCtrToPresent];
            end
            
            imgToPresent = find(ctrArray_Match==ctrMatchToPresent(end));
            imgMtrUpdated = bgc + (bgc*ctrArray_Match(imgToPresent)).*imgMtrTest{randi(nDynmAdpt)};
            imgTex_Match = Screen('MakeTexture',winPtr,imgMtrUpdated);
            Screen('DrawTexture',winPtr,imgTex_Match,[],imgRect_Match);
            
            imgMtrUpdated = bgc + (bgc*ctrArray_Test).*imgMtrTest{randi(nDynmAdpt)};
            imgTex_Test = Screen('MakeTexture',winPtr,imgMtrUpdated);            
            Screen('DrawTexture',winPtr,imgTex_Test,[],imgRect_Test);
            
            Screen('FillOval',winPtr,white,fixRect_w);
            Screen('FillOval',winPtr,black,fixRect_b);
            
            % control duration
            flipStartTime = Screen('Flip',winPtr);
            Rec_temp = [Rec_temp;testStartCode,ctrMatchToPresent(end),flipStartTime];
            while 1
                timeTmp = GetSecs;
                if (timeTmp - flipStartTime) > Dur_Test - flipAheadTime
                    break;
                end
                [keyIsDown,keySecs,keyCode] = KbCheck(-3);
                
                if keyIsDown == 0
                    isKeyRecordTest = 1;
                    isKeyRecordMatch = 1;
                end
                
                if keyIsDown && keyCode(keyTest) && isKeyRecordTest
                    Rec_temp = [Rec_temp;keyPrsCodeTest,0,keySecs];
                    isKeyRecordTest = 0;
                end                
                if keyIsDown && keyCode(keyMatch) && isKeyRecordMatch
                    Rec_temp = [Rec_temp;keyPrsCodeMatch,0,keySecs];
                    isKeyRecordMatch = 0;
                end                
                if keyIsDown && keyCode(keyEsc)
                    Rec = [Rec; Rec_temp];
                    Rec(:,end) = Rec(:,end) - stiStart;                                        
                    temp = clock; save(['ctrKeyRec_Cond',expCond,'_',sbjName,'_',num2str(temp(1:end-1))],'Rec','randseed');                    
                    RestoreCluts; Screen('CloseAll'); ShowCursor;ListenChar(0);
                    return;
                end
            end
            
            %% present Blank2
            Screen('FillOval',winPtr,white,fixRect_w);
            Screen('FillOval',winPtr,black,fixRect_b);
            
            % control duration
            flipStartTime = Screen('Flip',winPtr);
            while 1
                timeTmp = GetSecs;
                if (timeTmp - flipStartTime) > Dur_Blank2 - flipAheadTime
                    break;
                end
                [keyIsDown,keySecs,keyCode] = KbCheck(-3);
                
                if keyIsDown == 0
                    isKeyRecordTest = 1;
                    isKeyRecordMatch = 1;
                end
                                
                if keyIsDown && keyCode(keyTest) && isKeyRecordTest
                    Rec_temp = [Rec_temp;keyPrsCodeTest,0,keySecs];
                    isKeyRecordTest = 0;
                end
                if keyIsDown && keyCode(keyMatch) && isKeyRecordMatch
                    Rec_temp = [Rec_temp;keyPrsCodeMatch,0,keySecs];
                    isKeyRecordMatch = 0;
                end
                if keyIsDown && keyCode(keyEsc)
                    Rec = [Rec; Rec_temp];
                    Rec(:,end) = Rec(:,end) - stiStart;                                    
                    temp = clock; save(['ctrKeyRec_Cond',expCond,'_',sbjName,'_',num2str(temp(1:end-1))],'Rec','randseed');                    
                    RestoreCluts; Screen('CloseAll'); ShowCursor;ListenChar(0);
                    return;
                end
            end            
        end
        Rec = [Rec; Rec_temp];
        if expStage == '1'
            DrawFormattedText(winPtr,['Now you will go through 25 minutes of adaptation without rest.'],'center',winCenterY-30,white);
            DrawFormattedText(winPtr,['Please fixate on the center point but pay attention to the adaptors.'],'center',winCenterY+0,white);
            DrawFormattedText(winPtr,['Judge the contrast based on overall saliency.'],'center',winCenterY+30,white);
            DrawFormattedText(winPtr,['Press SPACE to start when you are ready.'],'center',winCenterY+60,white);
            flipStartTime = Screen('Flip',winPtr);
            while 1
                [keyIsDown,~,keyCode] = KbCheck(-3);
                if keyIsDown && keyCode(keySpace)
                    break
                end
            end
        end
        % a 1s blank after space key
        tempFlipStartTime = Screen('Flip',winPtr);
        Screen('Flip',winPtr,tempFlipStartTime+1);
        
    end
    Rec(:,end) = Rec(:,end) - stiStart;        
    temp = clock; save(['ctrKeyRec_Cond',expCond,'_',sbjName,'_',num2str(temp(1:end-1))],'Rec','randseed');
end

%% close window
RestoreCluts; Screen('CloseAll'); ShowCursor; ListenChar(0);
