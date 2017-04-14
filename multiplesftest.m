clear all; close all;
commandwindow;
smallwin = input('Small Window 1/0: ');
%% Global Var
adaptorAng = 20;
numPhases = 60;
siz = 600;
sf1 = 30;  %Cycles per stimulus
sf2 = 15;
gabWidth = 90;
frameDur = 0.1;
fixpix = 16; %size of the fixation cross
movieDurationSecs = 10;
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
%% Open Screen
if smallwin == 1
    [w,allcords] = Screen('OpenWindow',screenNumber, gray,[0,0,400,400]);
else
    [w,allcords] = Screen('OpenWindow',screenNumber, gray);
end
[xCenter,yCenter] = RectCenter(allcords); %get the center coordinates
fixRect1 = SetRect(0,0,fixpix,fixpix);
fixRect1 = CenterRect(fixRect1,allcords);
fixRect2 = SetRect(0,0,fixpix./2,fixpix./2);
fixRect2 = CenterRect(fixRect2,allcords);
%% Adaptor parameter
[x,y]=meshgrid(1:siz,1:siz);  %Makes an array of x and y values
x = x-siz./2; y=y-siz./2;
angle=adaptorAng*pi/180; %Convert degrees to radians
a=cos(angle);
b=sin(angle);
for i=1:numPhases
    phase=(i/numPhases)*2*pi;
    % Gabor
    gaus = exp(-((x/gabWidth).^2)-((y/gabWidth).^2));
    msf1=gaus.*sin(sf1*2*pi/siz*(a*x+b*y)+phase);
    msf2=gaus.*sin(sf2*2*pi/siz*(a*x+b*y)+phase);
    tex(i)=Screen('MakeTexture', w, round(gray+inc*msf1));
    tex2(i)=Screen('MakeTexture', w, round(gray+inc*msf2));
end
alltex = [tex tex2];
%% Movie Frames
movieDurationFrames=round(movieDurationSecs./frameDur);
% movieFrameIndices = randi(numPhases,movieDurationFrames,1);
movieFrameIndices = randi(numel(alltex),movieDurationFrames,1);

for i=1:movieDurationFrames
    Screen('DrawTexture', w, alltex(movieFrameIndices(i)));
%     Screen('DrawTexture', w, tex(movieFrameIndices(i)));
%     Screen('DrawTexture', w, tex2(movieFrameIndices(i))); %issue:tex is covered by tex2 
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

Screen('CloseAll');