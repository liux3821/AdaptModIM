clear all; close all;
smallwin = input('Small Window 1/0: ');

%% Global Var
testgrat1 = 45;
testgrat2 = 135;
siz = 600;
sf = 15;
gabwidth = 90;
fixpix = 16;
%% Call and Open Screen
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
[w,allcords] = Screen('OpenWindow',screenNumber, gray,[0,0,400,400]);
[xCenter,yCenter] = RectCenter(allcords); %get the center coordinates
fixRect1 = SetRect(0,0,fixpix,fixpix);
fixRect1 = CenterRect(fixRect1,allcords);
fixRect2 = SetRect(0,0,fixpix./2,fixpix./2);
fixRect2 = CenterRect(fixRect2,allcords);
%% Plaid parameter
[x,y]=meshgrid(1:siz,1:siz);  %Makes an array of x and y values
x = x-siz./2; y=y-siz./2;
angle1=testgrat1*pi/180; %Convert degrees to radians
angle2=testgrat2*pi/180;
a1=cos(angle1);
b1=sin(angle1);
a2=cos(angle2);
b2=sin(angle2);
%% Make Plaid
gaus = exp(-((x/gabwidth).^2)-((y/gabwidth).^2));
msf1=0.5.*(gaus.*sin(sf*2*pi/siz*(a1*x+b1*y)));
msf2=0.5.*(gaus.*sin(sf*2*pi/siz*(a2*x+b2*y))); %contrast*
allgrat = msf1+msf2;
allgrattex = Screen('MakeTexture',w,round(gray+inc*allgrat));
%tex45=Screen('MakeTexture', w, round(gray+inc*msf1));
%tex135=Screen('MakeTexture', w, round(gray+inc*msf2));
%Screen('DrawTexture', w, tex45);
%Screen('Flip', w);
%Screen('DrawTexture', w, tex135); 
%Screen('Flip', w);
Screen('DrawTexture', w, allgrattex);
% Screen('Flip', w);
Screen('FillOval', w, gray, fixRect1);
Screen('FillOval', w, black, fixRect2);
Screen('Flip', w);
WaitSecs(5);

Screen('CloseAll');