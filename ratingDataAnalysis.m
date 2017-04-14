
%CAT1:NATURAL SCENE
%CAT2:URBAN SCENE
%CAT3:BODY PARTS
%CAT4:TEXTURES
%CAT5:MANMADE OBJECTS
%CAT6:FACES

clear all; close all;
matpath = pwd;
commandwindow; 
sbjCode = input('File Name: ', 's');

resultFilePath = strcat(pwd,'/imgRatingResult/',sbjCode,'.mat');
load(resultFilePath);

tempOrder = []; ratingLabel = []; orderRating = [];
numofimages = length(imgIdx);

[sorted,sortIdx] = sort(imgIdx);
orderRating = imgRating(sortIdx);
stem(orderRating,'LineStyle','none') %axis([xmin xmax ymin ymax])
set(gca, 'YLim', [0,4])
line([30 30],[0 40],'Color',[0 0 0])
line([55 55],[0 40],'Color',[0 0 0])
line([78 78],[0 40],'Color',[0 0 0])
%line([115 115],[0 40],'Color',[0 0 0])
%line([145 145],[0 40],'Color',[0 0 0])