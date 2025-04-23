function motionA = maxxcorrAx(volume)

% volSize = size(volume); numFrame = volSize(3);
% 
% romA = 20;
% maskT=romA+1;
% maskB=volSize(1)-romA;
% motionA = zeros(1,numFrame);
% 
% for frameNum = 2:numFrame
%     mask = volume(maskT:maskB, 1:volSize(2), frameNum-1);
%     field = volume(:,:,frameNum);
%      for i = -romA:romA
%          fieldIn = field(maskT+i:maskB+i, :);   
%          cc(maskT+i,:) = corr2(mask, fieldIn);
%      end
%     [max_cc, imax] = max(abs(cc(:)));
%     [ypeak, xpeak] = ind2sub(size(cc),imax(1));
%     cc_offset = ypeak-romA-1;
%     motionA(frameNum) = cc_offset(1) + motionA(frameNum-1);
% end

volSize = size(volume); numFrame = volSize(3);

romA = 10;
maskT=romA+1;
maskB=volSize(1)-romA;
motionA = zeros(1,numFrame);

for i = 1:size(volume,3)
    [~, Gy] = imgradientxy((abs(volume(:,:,i)).^2));
    gradVolume(:,:,i) = Gy;
end

for frameNum = 2:numFrame
    mask = gradVolume(maskT:maskB, 1:volSize(2), frameNum-1);
    field = gradVolume(:,:,frameNum);
     for i = -romA:romA
         fieldIn = field(maskT+i:maskB+i, :);   
         cc(maskT+i,:) = corr2(mask, fieldIn);
     end
    [max_cc, imax] = max(abs(cc(:)));
    [ypeak, xpeak] = ind2sub(size(cc),imax(1));
    cc_offset = ypeak-romA-1;
    motionA(frameNum) = cc_offset(1) + motionA(frameNum-1);
end

end

