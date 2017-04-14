clear all; close all;
commandwindow;
smallwin = input('Small Window y/n: ', 's');
%% Global Variable
gratAngTask = 20;
gratAngRef = 0;
task_size = 4; %in degrees
reference_size = 4; %in degrees
cpd = 4;
cps_task = cpd*task_size;
cps_reference = cpd*reference_size;%Cycles per degree
gabWidth = 5;
fixpix = 8; %size of the fixation cross
maskSiz = 0.8;
num_trials = 5;
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
dpp = 180./(monitor.center(1).*2);
tasksizepixel = visAng2xyNew(task_size,0,monitor);
refsizepixel = visAng2xyNew(reference_size,0,monitor);

offsetpixleft = -1.*tasksizepixel; %85;
white=WhiteIndex(screenNumber);
black=BlackIndex(screenNumber);

gray=round((white+black)/2);

if gray == white
    gray=white / 2;
end

% Contrast 'inc'rement range for given white and gray values:
inc=white-gray;
newinc = 0.5.*inc;

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
imgRect_Adpt_Mid = CenterRectOnPoint([0,0,tasksizepixel,tasksizepixel],xCenter,yCenter);
imgRect_Task = CenterRectOnPoint([0,0,refsizepixel,refsizepixel],xCenter,yCenter);
imgRect_Reference = OffsetRect(imgRect_Adpt_Mid,offsetpixleft,0);

%TASK GRATING
taskGrating = drawGrating(gratAngTask,tasksizepixel,cps_task,0,gabWidth...
    ,maskSiz,0);
task_tex = Screen('MakeTexture', w, round(gray+inc*taskGrating));
%REFERENCE GRATING
refGrating = drawGrating(gratAngRef,refsizepixel,cps_reference,0,gabWidth...
    ,maskSiz,0);
ref_tex = Screen('MakeTexture', w, round(gray+newinc*refGrating));



%% Start Task
%Get Mouse Movement to adjust Ort
SetMouse(xCenter,yCenter);
[xMouse,yMouse,buttons] = GetMouse(w);

for iteration = 1:num_trials
    breakLoop = 0;
    while breakLoop == 0
        [adj_x,adj_y,buttons] = GetMouse(w);
        if adj_x > xCenter
            ort_change = ((adj_x - xCenter)*dpp)+gratAngTask;
        else
            ort_change = gratAngTask-((xCenter - adj_x)*dpp);
        end
        
        task_angle_change=ort_change*pi/180; %Convert degrees to radians
        task_grating_change = drawGrating(ort_change,tasksizepixel,...
            cps_task,0,gabWidth,maskSiz,0);  
        task_tex_change = Screen('MakeTexture', w, round(gray+inc*task_grating_change));
        
        
        timeStamp = GetSecs;
        while GetSecs - timeStamp <= 0.5 && breakLoop == 0
            Screen('DrawTexture', w, task_tex_change,[],imgRect_Task);
            Screen('DrawTexture', w, ref_tex,[],imgRect_Reference);
            Screen('FillOval', w, gray, fixRect1);
            Screen('FillOval', w, black, fixRect2);
            Screen('Flip',w);
            [xMouse,yMouse,buttons] = GetMouse(w);
            if buttons(1) == 1
                %save('result.mat','ort_change','-append');
                SetMouse(xCenter,yCenter);
                breakLoop = 1;
            end
        end
        
        timeStamp2 = GetSecs;
        while GetSecs - timeStamp2 <= 1 && breakLoop == 0
            grayTex = Screen('MakeTexture',w,gray);
            Screen('DrawTexture',w,grayTex);
            Screen('Flip',w);
            [xMouse,yMouse,buttons] = GetMouse(w);
            if buttons(1) == 1
                %save('result.mat','ort_change','-append');
                SetMouse(xCenter,yCenter);
                breakLoop = 1;
            end
        end
    end 
    Screen('DrawTexture',w,grayTex);
    Screen('Flip',w);
    WaitSecs(0.25);
end
%% Close Screens
Screen('Close');

RestoreCluts;
Screen('CloseAll'); ShowCursor; ListenChar(0);