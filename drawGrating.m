function [grat1,grat2] = drawGrating(ort, gratPixel, cps1,cps2,gabWidth, maskSiz,phase)

%ORT in degrees
%SF in cycles/deg
%Phase needed to define numPhases and put into a loop
%Otherwise set Phase = 0


maskLenGrating = maskSiz*gratPixel;

[GratingX,GratingY]=meshgrid(1:gratPixel,1:gratPixel);  %Makes an array of x and y values
GratingX = GratingX-gratPixel./2; GratingY=GratingY-gratPixel./2;

rGrating = sqrt(GratingX.^2+GratingY.^2);
maskGrating = rGrating<(maskLenGrating./2);

gratAngle=ort*pi/180;
angle_1=cos(gratAngle);
angle_2=sin(gratAngle);
gaus = exp(-((GratingX/gabWidth).^2)-((GratingY/gabWidth).^2));
maskGrating = conv2(gaus,double(maskGrating),'same');
maskGrating = maskGrating./max(max(maskGrating)); 
grat1=maskGrating.*sin(cps1*2*pi*(angle_1*GratingX+angle_2*GratingY)/gratPixel+phase);
grat2=maskGrating.*sin(cps2*2*pi*(angle_1*GratingX+angle_2*GratingY)/gratPixel+phase);

