clearvars
addpath(genpath('/Users/YM/Documents/GitHub/COIL/Project-PD_OCT/Melanoma/code'));
cd('C:\Users\coil_\Desktop\Github\Project-PD_OCT\Melanoma\code\JonesMatrixCode\JonesMatrixCode_Des_Ver6_ProcForPSSAO-OCT')
%%
figure(1), 
for i = 1:10:1000
    imshow(flip(1-DOPU_test(:,:,i)),[0 0.4]),colormap('gray')    
    pause(1)
end

entropy = flipud(1-DOPU_test);
% entropy(entropy>0)=1;

figure(2), 
for i = 1:10:1000
    imshow(entropy(:,:,i)),colormap('gray')    
    pause(1)
end

[DOPU_depth_C,DOPU_depth_C_bot]= genDOPU_depth_C(flipud(DOPU_test));
K = (1/10000)*ones(100);
DOPU_depth_C_sm = conv2(DOPU_depth_C,K,'same');
DOPU_depth_C_bot_sm = conv2(DOPU_depth_C_bot,K,'same');

figure(3), 
for i = 1:10:1000
    imshow(flip(1-DOPU_test(:,:,i)),[0 0.4]),colormap('hot'), hold on
    plot(DOPU_depth_C_sm(:,i)),plot(DOPU_depth_C_bot_sm(:,i)), hold off
    pause(0.1)
end

avgOCT = flipud((avgOCT));
% avgOCT = max(avgOCT(:))-avgOCT;

rHite = 5;
rDOPU = zeros(size(DOPU_depth_C));
rDOPU_vol = entropy;
for i = 1:1000
    for j = 1:1000
        upperlim = round(DOPU_depth_C_sm(i,j));
        lowerlim = upperlim + rHite;
        lowerlim(lowerlim>size(entropy,1)) = 1;
        upperlim(lowerlim>size(entropy,1)) = 1;
        rDOPU(i,j) = 1-mean(entropy(upperlim:lowerlim,i,j),1);
%         rDOPU(i,j) = 1-mean(avgOCT(upperlim:lowerlim,i,j),1);  
        rDOPU_vol(1:upperlim-1,i,j) = 0;
        rDOPU_vol(lowerlim+1:end,i,j) = 0;
    end
    if mod(i,100) == 0
        i
    end
end
figure,imshow(1-rDOPU,[]),colormap('hot')

savepath = '/Users/YM/Desktop/';
exportTiff(rDOPU_vol, fullfile(savepath,['rDOPU']))  