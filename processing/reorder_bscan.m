function vol_reorder = reorder_bscan(vol,numMscans)

% flip bscan %
vol_reorder = zeros(size(vol));
for FrameNum = 1:size(vol,3)
    if mod(FrameNum,numMscans*2) == 0 || mod(FrameNum,numMscans*2) > numMscans
        vol_reorder(:,:,FrameNum) = fliplr(vol(:,:,FrameNum));
    else
        vol_reorder(:,:,FrameNum) = vol(:,:,FrameNum);
    end
end