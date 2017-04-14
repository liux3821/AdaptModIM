function imgMtr = makeGrating(imgLen,ort,cycl,phs,maskLen,gsStr,isShow)
[x,y] = meshgrid(0:imgLen, 0:imgLen);
x = x-imgLen/2; y = y-imgLen/2;  % center at 0
x_theta=x*sin(ort) + y*cos(ort);
imgMtr = cos(2*pi/(imgLen/cycl)*x_theta + phs); % defines gratings [-1,1]

r = sqrt(x.^2 + y.^2);
imgMtr(r>(imgLen./2)) = 0; % defines circular aperture

% large gaussian variance: distributed
% small gaussian variance: centered
% large maskLen: goes off smoothly
% small maskLen: goes off quickly
gauss = normpdf(r,0,gsStr*imgLen/2);
mask = r<(maskLen./2);
mask = conv2(gauss,double(mask),'same');
mask = mask./max(max(mask));  % normalize mask

imgMtr = imgMtr.*mask; % convolution with Gaussian mask

if isShow
    figure; imagesc(imgMtr); colormap('gray');
end

end