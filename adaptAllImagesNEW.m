

%% HouseKeeping
clear all; close all;
KbName('UnifyKeyNames');
commandwindow; 
smallwin = input('Small Window y/n: ', 's');
testcond = input('Adaptor Same Orientation y/n: ', 's');
subject_code = input('Subject Code: ', 's');
exp_Loc = input('p for PsyPhys, l for Laptop: ', 's');
ListenChar(2); HideCursor;
alldir = pwd;
categories = {'urban_scenes','faces'};
allimg = [];
imgData = [];
imgRMS = [];
imgName = [];
imgRating = [];
imgIdx = [];
postImg = [];
tempimgName = []; wantSize = [450 450];

%% Adaptor Variables
% ORIENTATION
RMSCtrst = 0.3;
imgPresentDur = 0.3;
intLevDur = 0.2;
gratAngTest = [20, -20]; %EXPERIMENTAL CONDITION
gratAngCtrl = [20,20]; %CONTROL CONDITION
numPhases = 4;

%ANGULAR SIZE
%ATTN: ON LAPTOP SHOULD USE SMALLER ANGULAR SIZE:4.
%IN PSYCHOPHYSICS ROOM SHOULD USE ANGULAR SIZE 8.

adaptorsize = 8; %in degrees
testsize = 8; %in degrees
%imagesize = 500;

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
testgabWidth = 125; %IMAGE GABOR SIZE
frameDur = 0.1;
testDur = 1; %IMAGE PRESENT TIME
maskSiz = 0.8;
maskSizIMG = 0.8;
mask = maskSizIMG.*wantSize(1);

%EXPERIMENT PARAMETERS
adapttime = 2; %specify preadaptation duration
topUpBTTrial = 1;
fixpix = 8; %size of the fixation cross

%VECTORS STORE RESULTS
imageName = [];
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
    offsetpix = 105;
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

adaptorsizepixel = visAng2xyNew(adaptorsize,0,monitor);

%offsets sizes
offsetpixup = adaptorsizepixel-offsetpix; %95 pixels if on laptap
offsetpixdown = -1.*offsetpixup;

white=WhiteIndex(screenNumber);
black=BlackIndex(screenNumber);
gray=round((white+black)/2);

if gray == white
    gray=white / 2;
end

inc=white-gray; %100% MICHELSON CONTRAST
adaptCtrst = 0.9.*inc;%ADAPTOR CTRST:90%

%% Calibration
if exp_Loc == 'p'
    addpath([pwd,'/loadCalibration']);
    load calibration-0-13-Jan-2017.mat
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
imgRect_Img_Mid = CenterRectOnPoint([0,0,wantSize(1),wantSize(2)],xCenter,yCenter);
%imgRect_Test = CenterRectOnPoint([0,0,testsizepixel,testsizepixel],xCenter,yCenter);
imgRect_Adpt(1,:) = OffsetRect(imgRect_Adpt_Mid,0,offsetpixup);
imgRect_Adpt(2,:) = OffsetRect(imgRect_Adpt_Mid,0,offsetpixdown);
%DETERMINE GRATING ORIENTATIONS
if testcond == 'y'
    gratAng = gratAngCtrl;
else
    gratAng = gratAngTest;
end

%DRAW GRATING TEXTURE
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

grayTex = Screen('MakeTexture',w,gray);
% % Use realtime priority for better timing precision:
priorityLevel=MaxPriority(w);
Priority(priorityLevel);

%% READ ALL IMAGE DATA
for currdir=1:length(categories)
    tempimg = dir([alldir filesep 'testimages' filesep categories{currdir} filesep 'image*']);
    allimg = [allimg {tempimg.name}];
    numofimages = length(allimg);
end

for currimg = 1:length(allimg)
    tempim = allimg{currimg};
    whereformat = strfind(tempim,'.');
    catnum = str2num(tempim(whereformat-1));
    tempdata = double(imread([alldir filesep 'testimages' filesep categories{catnum} filesep allimg{currimg}]));
    %GREYSCALE IMAGE
    if (length(size(tempdata)) > 2)
        %NormImage = rgb2gray(NormImage);
        tempdata = mean(tempdata,3);
    else
        tempdata = tempdata;
    end
    imgData{currimg} = tempdata;
end
%RANDOMLLY DRAW IMAGE W/O REPLAECMENT
[randomImgCell,imgIdx] = datasample(imgData,numofimages,'Replace',false);

for counter = 1:length(imgIdx)
    tempimgName = allimg{imgIdx(counter)};
    imgName{counter} = tempimgName;
end

%IMAGE BLURR
[imgX,imgY] = meshgrid(1:wantSize(2), 1:wantSize(1));
%[kernX,kernY] = meshgrid(1:50, 1:50);
%kernX = kernX-50./2; kernY = kernY-50./2;
imgX = imgX-wantSize(2)./2; imgY = imgY-wantSize(1)./2;
%rImg = sqrt(imgX.^2+imgY.^2);
%maskIMG = rImg<(mask./2);
imageGaus = exp(-((imgX/testgabWidth).^2)-((imgY/testgabWidth).^2));
%imageGaus = exp(-((kernX/testgabWidth).^2)-((kernY/testgabWidth).^2));
% maskImage = conv2(imageGaus, double(maskIMG), 'same');
%maskImage = conv2(double(maskIMG), imageGaus, 'same');
%maskImage = maskImage./max(max(maskImage));
for i = 1:numofimages
    randomImgMat = randomImgCell{i};
    [NormedImg] = normImageNEW(randomImgMat,RMSCtrst,wantSize,testgabWidth);
    NormedImg(:,:,2) = 255.*imageGaus; %255 is 0 transparency
    tempNormIMG = NormedImg;
    postImg{i} = tempNormIMG;
end

%% PROMPT SCREEN
StartDEMO = 0;
while ~StartDEMO
    DrawFormattedText(w,['\n\n\n\n\n\n\n\n'...
        'Press RETURN to HAVE SOME FUN'],'center','center',white);
    Screen('FillOval', w, gray, fixRect1);
    Screen('FillOval', w, black, fixRect2);
    Screen('Flip',w);
    [KeyDown,time,keyCode] = KbCheck;
    if KeyDown == 1
        if keyCode(40)
            StartDEMO = 1;
            break;
        end
    end
end

%% START ADAPTATION
%adptStartTimeStamp = GetSecs; 
%trial_num = 1;
for iteration = 1:numofimages
    if iteration == 1
        movieDurationSecs = adapttime;
    else
        movieDurationSecs = topUpBTTrial;
    end
    
    movieDurationFrames=round(movieDurationSecs./frameDur);
    movieFrameIndices = randi(numPhases,movieDurationFrames,1);
    
    for i=1:movieDurationFrames
        for updown = 1:2
            Screen('DrawTexture', w, tex(updown, movieFrameIndices(i)),[],imgRect_Adpt(updown,:));
        end    
        Screen('FillOval', w, gray, fixRect1);
        Screen('FillOval', w, black, fixRect2);
        Screen('Flip', w);
        adaptTimeStamp = GetSecs;
        while GetSecs - adaptTimeStamp <= 0.1
            if iteration > 1
                [KeyDown,time,keyCode] = KbCheck(-3);
                if KeyDown == 1
                    if keyCode(80)
                        tempRating = 1;
                        imgRating(iteration-1) = tempRating;
                    elseif keyCode(81)
                        tempRating = 2;
                        imgRating(iteration-1) = tempRating;
                    elseif keyCode(79)
                        tempRating = 3;
                        imgRating(iteration-1) = tempRating;
                    else
                        imgRating(iteration-1) = nan;
                    end
                    
                end
            end
        end

    end
    if length(imgRating) < (iteration-1) && iteration > 1
        imgRating(iteration-1) = 998;
    end
    
    blankStart = GetSecs;
    Screen('DrawTexture',w,grayTex);
    while GetSecs - blankStart <= intLevDur
        Screen('FillOval', w, gray, fixRect1);
        Screen('FillOval', w, black, fixRect2);
        Screen('Flip',w);
    end
    
    %randomImgMat = randomImgCell{iteration};
    %[NormedImg] = normImageNEW(randomImgMat,RMSCtrst,wantSize,testgabWidth);
    %NormedImg(:,:,2) = 255.*imageGaus; %255 is 0 transparency
    %imgIT = num2str(iteration);
    %{
    DrawFormattedText(w,['\n\n\n\n\n\n\n\n'...
        imgIT],'center','center',white);
    %}
    testtex = Screen('MakeTexture', w, postImg{iteration});
    imgStart = GetSecs;
    while GetSecs - imgStart <= imgPresentDur
        Screen('DrawTexture', w, testtex,[],imgRect_Img_Mid);
        Screen('FillOval', w, gray, fixRect1);
        Screen('FillOval', w, black, fixRect2);
        Screen('Flip', w);
    end
    
    Screen('Close',testtex);
    
    blankStart2 = GetSecs;
    while GetSecs - blankStart2 <= intLevDur
        Screen('DrawTexture',w,grayTex);
        Screen('FillOval', w, gray, fixRect1);
        Screen('FillOval', w, black, fixRect2);
        Screen('Flip',w);
    end
end

%% END RESULT
for i=1:movieDurationFrames
    
    movieDurationFrames=round(topUpBTTrial./frameDur);
    movieFrameIndices = randi(numPhases,movieDurationFrames,1);
    
    for updown = 1:2
        Screen('DrawTexture', w, tex(updown, movieFrameIndices(i)),[],imgRect_Adpt(updown,:));
    end
    Screen('FillOval', w, gray, fixRect1);
    Screen('FillOval', w, black, fixRect2);
    Screen('Flip', w);
    adaptTimeStamp = GetSecs;
    while GetSecs - adaptTimeStamp <= 0.1
        [KeyDown,time,keyCode] = KbCheck(-3);
        if KeyDown == 1
            if keyCode(80)
                tempRating = 1;
                imgRating(numofimages) = tempRating;
            elseif keyCode(81)
                tempRating = 2;
                imgRating(numofimages) = tempRating;
            elseif keyCode(79)
                tempRating = 3;
                imgRating(numofimages) = tempRating;
            else
                imgRating(numofimages) = nan;
            end
            
        end
    end
end

%% Save Test Result
resultFilePath = strcat(pwd,'/imgRatingResult/',subject_code,'.mat');
clear imgData; clear tempdata; clear randomImgCell; clear randomImgMat; clear NormedImg; clear imageGaus;
save(resultFilePath);

%% Close Screens
Screen('Close');
RestoreCluts;
Screen('CloseAll'); ShowCursor; ListenChar(0);
