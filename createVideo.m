close all
clear all


%%

% 
rootImgFolder = 'Images';
outFolder = 'SimulationVideos';
outName = ['solution_narrowed',datestr(datetime('now'))];

fileStruct = dir(rootImgFolder);
fileIndexes = 3:length(fileStruct);

duration = 6;

v = VideoWriter([outFolder,'/',outName,'.avi'],'Motion JPEG AVI');
v.Quality = 85;
v.FrameRate = ceil(length(fileIndexes)/duration);
open(v);


for processIdx = 1:(length(fileStruct)-2);


    imgName = fileStruct(fileIndexes(processIdx)).name;
   
    M = imread( [rootImgFolder,'/',imgName]);

    writeVideo(v,M);
    
    
end

close(v);