%function XinyuDemo
% DriftDemo

gratAng = 15;
numPhases = 10;
siz = 600;
f = 15;  %Cycles per stimulus
gabWidth = 90;
frameDur = 0.1;
movieDurationSecs=50;

img = imread('test.png');
img = imresize(img, [100 300]);

try
	% This script calls Psychtoolbox commands available only in OpenGL-based 
	% versions of the Psychtoolbox. (So far, the OS X Psychtoolbox is the
	% only OpenGL-base Psychtoolbox.)  The Psychtoolbox command AssertPsychOpenGL will issue
	% an error message if someone tries to execute this script on a computer without
	% an OpenGL Psychtoolbox
	AssertOpenGL;
	
	% Get the list of screens and choose the one with the highest screen number.
	% Screen 0 is, by definition, the display with the menu bar. Often when 
	% two monitors are connected the one without the menu bar is used as 
	% the stimulus display.  Chosing the display with the highest dislay number is 
	% a best guess about where you want the stimulus displayed.  
	screens=Screen('Screens');
	screenNumber=max(screens);
	
    % Find the color values which correspond to white and black: Usually
	% black is always 0 and white 255, but this rule is not true if one of
	% the high precision framebuffer modes is enabled via the
	% PsychImaging() commmand, so we query the true values via the
	% functions WhiteIndex and BlackIndex:
	white=WhiteIndex(screenNumber);
	black=BlackIndex(screenNumber);
    
    % Round gray to integral number, to avoid roundoff artifacts with some
    % graphics cards:
	gray=round((white+black)/2);

    % This makes sure that on floating point framebuffers we still get a
    % well defined gray. It isn't strictly neccessary in this demo:
    if gray == white
		gray=white / 2;
    end
    
    % Contrast 'inc'rement range for given white and gray values:
	inc=white-gray;

    % Open a double buffered fullscreen window and select a gray background 
	% color:
	w=Screen('OpenWindow',screenNumber, gray);
    
	% Compute each frame of the movie and convert the those frames, stored in
	% MATLAB matices, into Psychtoolbox OpenGL textures using 'MakeTexture';
    
    %Grating parameter
    	[x,y]=meshgrid(1:siz,1:siz);  %Makes an array of x and y values
        x = x-siz./2; y=y-siz./2;
		angle=gratAng*pi/180; %Convert degrees to radians
		a=cos(angle);
		b=sin(angle);

	for i=1:numPhases
		phase=(i/numPhases)*2*pi;
		% Gabor
        m=exp(-((x/gabWidth).^2)-((y/gabWidth).^2)).*sin(f*2*pi/siz*(a*x+b*y)+phase);
		tex(i)=Screen('MakeTexture', w, gray+inc*m); %#ok<AGROW>
    end
		
    testtex=Screen('MakeTexture', w, img);

	% Run the movie animation for a fixed period.  
    % Convert movieDuration in seconds to duration in frames to draw:
    movieDurationFrames=round(movieDurationSecs./frameDur);
    
    %Put random phase grating
	%movieFrameIndices=mod(0:(movieDurationFrames-1), numPhases) + 1;  %Play straight through
    movieFrameIndices = randi(numPhases,movieDurationFrames,1);
    
    % Use realtime priority for better timing precision:
    priorityLevel=MaxPriority(w);
	Priority(priorityLevel);

    % Animation loop:
    for i=1:movieDurationFrames
        % Draw image:
        Screen('DrawTexture', w, tex(movieFrameIndices(i)));
        % Show it at next display vertical retrace. Please check DriftDemo2
        % and later, as well as DriftWaitDemo for much better approaches to
        % guarantee a robust and constant animation display timing! This is
        % very basic and not best practice!
        Screen('Flip', w);
        WaitSecs(0.1);
    end

    Screen('DrawTexture', w, testtex);
    Screen('Flip', w);
    WaitSecs(2);
    
    Priority(0);
	
    % Close all textures. This is not strictly needed, as
    % Screen('CloseAll') would do it anyway. However, it avoids warnings by
    % Psychtoolbox about unclosed textures. The warnings trigger if more
    % than 10 textures are open at invocation of Screen('CloseAll') and we
    % have 12 textues here:
    Screen('Close');
    
    % Close window:
    Screen('CloseAll');

catch
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    Screen('CloseAll');
    Priority(0);
    psychrethrow(psychlasterror);
end %try..catch..
