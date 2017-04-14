clear all; close all;
commandwindow;
smallwin = input('Small Window 1/0: ');
%% Global Variable 
gratAng = -20;
testAng = 0;
numPhases = 10;
adaptorsize = 2; %in degrees
testsize = 8; %in degrees
cpd = 4;
cps_adaptor = cpd*adaptorsize; 
cps_test = cpd*testsize;%Cycles per degree
gabWidth = 20;
testWidth = 80;
frameDur = 0.1;
testDur = 1;
adapttime = 120; %specify preadaptation duration
topup = 8;
fixpix = 8; %size of the fixation cross
numoftrials = 10;
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

white=WhiteIndex(screenNumber);
black=BlackIndex(screenNumber);

gray=round((white+black)/2);

if gray == white
    gray=white / 2;
end

% Contrast 'inc'rement range for given white and gray values:
inc=white-gray;

%% Calibration
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
if smallwin == 1
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
imgRect_Adpt = CenterRectOnPoint([0,0,adaptorsizepixel,adaptorsizepixel],xCenter,yCenter);
imgRect_Test = CenterRectOnPoint([0,0,testsizepixel,testsizepixel],xCenter,yCenter);
%% Grating and test parameter
[adaptorx,adaptory]=meshgrid(1:adaptorsizepixel,1:adaptorsizepixel);  %Makes an array of x and y values
adaptorx = adaptorx-adaptorsizepixel./2; adaptory=adaptory-adaptorsizepixel./2;
angle=gratAng*pi/180; %Convert degrees to radians
adaptor1=cos(angle);
adaptor2=sin(angle);
%Test Parameter big
[testx,testy] = meshgrid(1:testsizepixel, 1:testsizepixel);
testx = testx-testsizepixel./2; testy = testy-testsizepixel./2;
testangle=testAng*pi/180;
test1=cos(testangle);
test2=sin(testangle);

for i=1:numPhases
    phase=(i/numPhases)*2*pi;
    % Gabor
    gaus = exp(-((adaptorx/gabWidth).^2)-((adaptory/gabWidth).^2));
    m=gaus.*sin(cps_adaptor*2*pi*(adaptor1*adaptorx+adaptor2*adaptory)/adaptorsizepixel+phase);
    tex(i)=Screen('MakeTexture', w, round(gray+inc*m)); %#ok<AGROW>
end

% Use realtime priority for better timing precision:
priorityLevel=MaxPriority(w);
Priority(priorityLevel);
%% Animation Loop
for iteration = 1:numoftrials
    
    if iteration == 1
        movieDurationSecs = adapttime;
    else
        movieDurationSecs = topup;
    end
    movieDurationFrames=round(movieDurationSecs./frameDur);
    movieFrameIndices = randi(numPhases,movieDurationFrames,1);
    for i=1:movieDurationFrames
        Screen('DrawTexture', w, tex(movieFrameIndices(i)),[],imgRect_Adpt);
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
 
    if mod(iteration,2) == 0
        testgaus = exp(-((testx/testWidth).^2)-((testy/testWidth).^2));
        m2=testgaus.*sin(cps_test*2*pi/testsizepixel*(test1*testx+test2*testy));
        gratingtex = Screen('MakeTexture', w, round(gray+inc*m2));
        Screen('DrawTexture', w, gratingtex,[],imgRect_Test);
        Screen('FillOval', w, gray, fixRect1);
        Screen('FillOval', w, black, fixRect2);
        Screen('Flip',w);
        WaitSecs(testDur);
    else
        mreg=gaus.*sin(cps_adaptor*2*pi/adaptorsizepixel*(test1*adaptorx+test2*adaptory));
        gratingtexreg = Screen('MakeTexture', w, round(gray+inc*mreg));
        Screen('DrawTexture', w, gratingtexreg,[],imgRect_Adpt);
        Screen('FillOval', w, gray, fixRect1);
        Screen('FillOval', w, black, fixRect2);
        Screen('Flip',w);
        WaitSecs(testDur);
    end
end

Screen('Close');

RestoreCluts;
Screen('CloseAll'); ShowCursor; ListenChar(0);
