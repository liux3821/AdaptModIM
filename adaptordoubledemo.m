clear all; close all;
commandwindow;
smallwin = input('Small Window y/n: ', 's');
testcond = input('Adaptor Same Orientation y/n: ', 's');
%% Global Variable
gratAngTest = [20, -20];
gratAngCtrl = [20,20];
testAng = 0;
numPhases = 4;
adaptorsize = 4; %in degrees
testsize = 10; %in degrees
cpd = 4;
cpd_random = 6;
cps_adaptor = cpd*adaptorsize;
cps_test = cpd*testsize;
cps_random = cpd_random*adaptorsize;%Cycles per degree
gabWidth = 5;
gabtestgabWidth = 100;
testgabWidth = 150;
frameDur = 0.1;
testDur = 1;
adapttime = 5; %specify preadaptation duration
topup = 3;
fixpix = 8; %size of the fixation cross
numoftrials = 10;
maskSiz = 0.8;
randomStartFrame = 2;
randomFinishFrame = 28;
randomFrame = randi([randomStartFrame, randomFinishFrame]);
%% Call Screen
AssertOpenGL;
Screen('Preference', 'VisualDebugLevel', 1); %get rid of the welcoming screen
screens=Screen('Screens');
screenNumber=max(screens);
screenRect=Screen('Rect', screenNumber);
monitor.viewDist = 58; %cm
monitor.center = screenRect(3:4)./2;
monitor.size(1) = 29; %Also cm
monitor.size(2)= monitor.size(1).*(monitor.center(2)./monitor.center(1));

adaptorsizepixel = visAng2xyNew(adaptorsize,0,monitor);
testsizepixel = visAng2xyNew(testsize,0,monitor);

offsetpixup = adaptorsizepixel-95; %85;
offsetpixdown = -1.*offsetpixup;
white=WhiteIndex(screenNumber);
black=BlackIndex(screenNumber);

gray=round((white+black)/2);

if gray == white
    gray=white / 2;
end

% Contrast 'inc'rement range for given white and gray values:
inc=white-gray;

%% Calibration
addpath([pwd,'/calibrationResults']);
load calibration-0-27-Oct-2016.mat;
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
imgRect_Test = CenterRectOnPoint([0,0,testsizepixel,testsizepixel],xCenter,yCenter);
imgRect_Adpt(1,:) = OffsetRect(imgRect_Adpt_Mid,0,offsetpixup);
imgRect_Adpt(2,:) = OffsetRect(imgRect_Adpt_Mid,0,offsetpixdown);
%imgRect_Test = OffsetRect(imgRect_Test,0,testsizepixel);


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
            inc*randomGrating));
        tex(updown,i)=Screen('MakeTexture', w, round(gray+...
            inc*adaptorGrating));
    end
end
d = dir('/Users/liux3821/Dropbox/MATLAB/Engel_Lab/Adaptation_TAE/testimages/image_vert*');
numofimages = length(d);
% % Use realtime priority for better timing precision:
priorityLevel=MaxPriority(w);
Priority(priorityLevel);
%% Animation Loop
for iteration = 1:numofimages+1
    if iteration == 1
        movieDurationSecs = adapttime;
    else
        movieDurationSecs = topup;
    end
    movieDurationFrames=round(movieDurationSecs./frameDur);
    movieFrameIndices = randi(numPhases,movieDurationFrames,1);
    for i=1:movieDurationFrames
        for updown = 1:2
            if i == randomFrame
                Screen('DrawTexture', w, tex(updown, movieFrameIndices(i)),[],imgRect_Adpt(updown,:));
                Screen('DrawTexture',w,RandomGratingtex(updown, movieFrameIndices(i))...
                    ,[],imgRect_Adpt(updown,:));
            else
                Screen('DrawTexture', w, tex(updown, movieFrameIndices(i)),[],imgRect_Adpt(updown,:));
            end
        end
        Screen('FillOval', w, gray, fixRect1);
        Screen('FillOval', w, black, fixRect2);
        Screen('Flip', w);
        WaitSecs(0.1);
        [keyIsDown,foo,keyCode] = KbCheck;
        if (keyIsDown == 1)
            if keyCode(KbName('ESCAPE')) %
                Priority(0);
                Screen('Close');
                RestoreCluts;
                Screen('CloseAll');
                return
            end
        end
    end
    if iteration == 1
        testGrating = drawGrating(testAng,testsizepixel,cps_test,...
            0,testgabWidth,maskSiz,0);
        testGratingTex = Screen('MakeTexture', w, round(gray+inc*testGrating));
        Screen('DrawTexture', w, testGratingTex,[],imgRect_Test);
        Screen('FillOval', w, gray, fixRect1);
        Screen('FillOval', w, black, fixRect2);
        Screen('Flip',w);
        WaitSecs(testDur);
    else
        imageFilename = ['/Users/liux3821/Dropbox/MATLAB/Engel_Lab/Adaptation_TAE/testimages/' d(iteration-1).name];
        imageData = imread(imageFilename);
        if (length(size(imageData)) > 2)
            imageGray = rgb2gray(imageData);
        else
            imageGray = imageData;
        end
        %resize images
        newsize = [testsizepixel,testsizepixel];
        [test_x,test_y] = meshgrid(1:testsizepixel,1:testsizepixel);
        test_x = test_x-testsizepixel./2; test_y = test_y-testsizepixel./2;
        testimages = imresize(imageGray,newsize);
        testgaus = exp(-((test_x/testgabWidth).^2)-((test_y/testgabWidth).^2));
        testimages(:,:,2) = 255.*testgaus; %255 is 0 transparency
        %make gabor envelope
        %meanimages = double(mean(testimages(:)));
        %testgaus2 = exp(-((testx/testWidth).^2)-((testy/testWidth).^2));
        %gaborimages=((testimages-128))./(128).*testgaus2.*128+128;
        %gaborimages = (testgaus.*testimages);
        testtex = Screen('MakeTexture', w, testimages);
        Screen('DrawTexture', w, testtex);
        Screen('FillOval', w, gray, fixRect1);
        Screen('FillOval', w, black, fixRect2);
        Screen('Flip', w);
        WaitSecs(testDur);
    end
end

Screen('Close');

RestoreCluts;
Screen('CloseAll'); ShowCursor; ListenChar(0);
