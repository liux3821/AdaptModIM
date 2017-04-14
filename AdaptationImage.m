%Demo TAE in natural images
clear all; close all;
commandwindow;
adapt_ori = input ('Adaptation Orientation 0 for vertical 1 for horizontal: ');
smallwin = input('Small Window 1/0: ');
%% Global Var
%CHANGE ADAPTOR ORIENTATION: line5,6,15,137.
gratAngVert = 25;
gratAngHori = 75;
testAngVert = 0;
testAngHori = 90;
numPhases = 10;
siz = 600;
testsize = 600;
f = 30;  %Cycles per stimulus
gabWidth = 90;
testWidth = 150;
frameDur = 0.1;
adapttime = 100; %specify preadaptation duration
fixpix = 16; %size of the fixation cross

%% Call Screen
AssertOpenGL;
Screen('Preference', 'VisualDebugLevel', 1); %get rid of the welcoming screen
screens=Screen('Screens');
screenNumber=max(screens);

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
fixRect1 = SetRect(0,0,fixpix,fixpix);
fixRect1 = CenterRect(fixRect1,allcords);
fixRect2 = SetRect(0,0,fixpix./2,fixpix./2);
fixRect2 = CenterRect(fixRect2,allcords);

%% Grating and test patch parameter
if adapt_ori == 0
    gratAng = gratAngVert;
    testAng = testAngVert;
    d = dir('/Users/liux3821/Documents/MATLAB/Engel_Lab/Adaptation_TAE/testimages/image_vert*');
    numofimages = length(d);
else
    gratAng = gratAngHori;
    testAng = testAngHori;
    d = dir('/Users/liux3821/Documents/MATLAB/Engel_Lab/Adaptation_TAE/testimages/image_hori*');
    numofimages = length(d);
end
%Grating parameter
[x,y]=meshgrid(1:siz,1:siz);  %Makes an array of x and y values
x = x-siz./2; y=y-siz./2;
angle=gratAng*pi/180; %Convert degrees to radians
a=cos(angle);
b=sin(angle);

%test patch parameter
[testx,testy] = meshgrid(1:testsize, 1:testsize);
testx = testx-testsize./2; testy = testy-testsize./2;
testangle=testAng*pi/180;
test1=cos(testangle);
test2=sin(testangle);

for i=1:numPhases
    phase=(i/numPhases)*2*pi;
    % Gabor
    gaus = exp(-((x/gabWidth).^2)-((y/gabWidth).^2));
    m=gaus.*sin(f*2*pi/siz*(a*x+b*y)+phase);
    tex(i)=Screen('MakeTexture', w, round(gray+inc*m)); %#ok<AGROW>
end

% Use realtime priority for better timing precision:
priorityLevel=MaxPriority(w);
Priority(priorityLevel);

%% Movie loop
for iteration = 1:numofimages
    % decide whether preadapt or topup adapt (topup adapt is 15s)
    if iteration == 1
        movieDurationSecs = adapttime;
    else
        movieDurationSecs = 1;
    end
    movieDurationFrames=round(movieDurationSecs./frameDur);
    movieFrameIndices = randi(numPhases,movieDurationFrames,1);
    %% Animation loop:
    for i=1:movieDurationFrames
        % Draw image:
        Screen('DrawTexture', w, tex(movieFrameIndices(i)));
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
    
    %{
    showtex = Screen('MakeTexture', w, round(gray+inc*m));
    Screen('DrawTexture', w, showtex);
    Screen('FillOval', w, gray, fixRect1);
    Screen('FillOval', w, black, fixRect2);
    Screen('Flip', w, showtex);
    WaitSecs(5);
    %}
    
    if iteration == 1
        testgaus = exp(-((testx/gabWidth).^2)-((testy/gabWidth).^2));
        m2=testgaus.*sin(f*2*pi/testsize*(test1*testx+test2*testy));
        gratingtex = Screen('MakeTexture', w, round(gray+inc*m2));
        Screen('DrawTexture', w, gratingtex);
        Screen('Flip',w);
        Screen('FillOval', w, gray, fixRect1);
        Screen('FillOval', w, black, fixRect2);
        WaitSecs(2);
    else
        %% readin images
        imageFilename = ['/Users/liux3821/Documents/MATLAB/Engel_Lab/Adaptation_TAE/testimages/' d(iteration-1).name];
        %imageFilename = ['/Users/liux3821/Documents/MATLAB/Engel_Lab/Adaptation_TAE/testimages_vertical/image' num2str(iteration-1) '.jpg'];
        imageData = imread(imageFilename);
        if (length(size(imageData)) > 2)
            imageGray = rgb2gray(imageData);
        else
            imageGray = imageData;
        end
        %resize images
        newsize = [testsize,testsize];
        testimages = imresize(imageGray,newsize);
        testimages(:,:,2) = 255.*testgaus; %255 is 0 transparency
        %make gabor envelope
        %meanimages = double(mean(testimages(:)));
        testgaus2 = exp(-((testx/testWidth).^2)-((testy/testWidth).^2));
        %gaborimages=((testimages-128))./(128).*testgaus2.*128+128;
        %gaborimages = (testgaus.*testimages);
        testtex = Screen('MakeTexture', w, testimages);
        Screen('DrawTexture', w, testtex);
        Screen('FillOval', w, gray, fixRect1);
        Screen('FillOval', w, black, fixRect2);
        Screen('Flip', w);
        WaitSecs(2);
    end
end

Priority(0);

Screen('Close');

RestoreCluts;
Screen('CloseAll'); ShowCursor; ListenChar(0);




