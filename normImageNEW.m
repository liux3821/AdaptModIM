function[NormalizedImage] = normImageNEW(OldImageStat,RMSCtrst,wantSize,testgabWidth)

blankIm = ones(wantSize).*128;
edgeMask = zeros(wantSize);
NormalizedImage = imresize(OldImageStat, [wantSize(1) NaN]);


if size(NormalizedImage,2) > wantSize(2)
    rect = CenterRectOnPoint([0 0 wantSize(2) wantSize(1)], round(size(NormalizedImage,2)./2), wantSize(1)./2);
    NormalizedImage = NormalizedImage(:,rect(1)+1:rect(3));
else
    rect = CenterRectOnPoint([0 0 size(NormalizedImage,2) size(NormalizedImage,1)], wantSize(2)./2, wantSize(1)./2);
    blankIm(:,rect(1)+1:rect(3)) = NormalizedImage;
    edgeMask(:,rect(1)+1:rect(3)) = 1;
    NormalizedImage = blankIm;
end

%{
%EdgeBlurr Image
tic
[imgX,imgY] = meshgrid(1:wantSize(2), 1:wantSize(1));
[kernX,kernY] = meshgrid(1:50, 1:50);
kernX = kernX-50./2; kernY = kernY-50./2;
imgX = imgX-wantSize(2)./2; imgY = imgY-wantSize(1)./2;
rImg = sqrt(imgX.^2+imgY.^2);
maskIMG = rImg<(mask./2);
% imageGaus = exp(-((imgX/testgabWidth).^2)-((imgY/testgabWidth).^2));
imageGaus = exp(-((kernX/testgabWidth).^2)-((kernY/testgabWidth).^2));
% maskImage = conv2(imageGaus, double(maskIMG), 'same');
maskImage = conv2(double(maskIMG), imageGaus, 'same');
maskImage = maskImage./max(max(maskImage));
toc
%}

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
    