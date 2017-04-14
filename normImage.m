function[NormalizedImage,imageGaus] = normImage(OldImageStat,imgsize,RMSCtrst,testgabWidth,wantSize)

%wantSize1 = [768 576];
%wantSize2 = [576 768];
blankIm = ones(wantSize).*128;
NormalizedImage = imresize(OldImageStat, [wantSize(1) NaN]);
if size(NormalizedImage,2) > wantSize(2)
    rect = CenterRectOnPoint([0 0 wantSize(2) wantSize(1)], round(size(NormalizedImage,2)./2), wantSize(1)./2);
    NormalizedImage = NormalizedImage(:,rect(1)+1:rect(3));
else
    rect = CenterRectOnPoint([0 0 size(NormalizedImage,2) size(NormalizedImage,1)], wantSize(2)./2, wantSize(1)./2);
    blankIm(:,rect(1)+1:rect(3)) = NormalizedImage;
    NormalizedImage = blankIm;
end

[imgX,imgY] = meshgrid(1:wantSize(2), 1:wantSize(1));
imgX = imgX-wantSize(2)./2; imgY = imgY-wantSize(1)./2;
imageGaus = exp(-((imgX/testgabWidth).^2)-((imgY/testgabWidth).^2));

%GREYSCALE IMAGE
if length(size(NormalizedImage)) > 2
% NormImage = rgb2gray(NormImage);
NormalizedImage = mean(NormalizedImage,3);
    else
        NormalizedImage = NormalizedImage;
end

%RMS CONTRAST
%RMS (IMAGE CONTRAST)/RMS IDEAL VALUE
%MULTIPLY BY IMAGE CONTRAST
%CONVERT BACK TO PIXELS: (CONTRAST*MEAN)+MEAN
meanIMGCtrst = mean(NormalizedImage(:));
IMGCtrst = (NormalizedImage - meanIMGCtrst)./meanIMGCtrst;
RMSImage = rms(IMGCtrst(:));
RMSRatio = RMSCtrst./RMSImage;
RMSNormed = IMGCtrst.*RMSRatio;
%NormalizedImage = (RMSNormed.* meanIMGCtrst)+meanIMGCtrst;
NormalizedImage = (RMSNormed.* 128)+128;


% NormImage = 128.*(RMSCtrst.*IMGCtrst./RMSImage)+128;
%NormalizedImage = 0.5.*((RMSCtrst.*IMGCtrst./RMSImage)+1);
    