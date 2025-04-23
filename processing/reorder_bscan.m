% function vol_reorder = reorder_bscan(vol,numMscans)
% 
% % flip bscan %
% vol_reorder = zeros(size(vol));
% for FrameNum = 1:size(vol,3)
%     if mod(FrameNum,numMscans*2) == 0 || mod(FrameNum,numMscans*2) > numMscans
%         vol_reorder(:,:,FrameNum) = fliplr(vol(:,:,FrameNum));
%     else
%         vol_reorder(:,:,FrameNum) = vol(:,:,FrameNum);
%     end
% end

function vol_reorder = reorder_bscan(vol,numMscans)

vol_reorder = vol;
if numMscans == 1
    for FrameNum = 1:size(vol,3)
        if mod(FrameNum,numMscans*2) == 0 || mod(FrameNum,numMscans*2) > numMscans
            vol_reorder(:,:,FrameNum) = fliplr(vol(:,:,FrameNum));
        else
            vol_reorder(:,:,FrameNum) = vol(:,:,FrameNum);
        end
    end
elseif numMscans == 4 % VISTA 2 FBFB FBFB FB... = 0101 0101 23...
    fwdScans = zeros(size(vol,1),size(vol,2),size(vol,3)/2);
    bwdScans = zeros(size(vol,1),size(vol,2),size(vol,3)/2);
     f_idx = 1;
     b_idx = 1;
    for FrameNum = 1:size(vol,3)
            if mod(FrameNum,2) == 0
                bwdScans(:,:,b_idx) = fliplr(vol(:,:,FrameNum)); % need to flip the backward scans
                b_idx = b_idx +  1;
            else
                fwdScans(:,:,f_idx) = vol(:,:,FrameNum);
                f_idx = f_idx +  1;
            end
    end
    
    j=1;
    k=1;
    for i=1:numMscans*2:size(vol,3)-4  
        vol_reorder(:,:,i:i+3) = fwdScans(:,:,j:j+3);
%         if(i<size(vol,3)-4)
            j=j+4;
%         end
    end
    for i=numMscans+1:numMscans*2:size(vol,3)-4
        vol_reorder(:,:,i:i+3) = bwdScans(:,:,k:k+3);
        k=k+4;
    end
end

%     elseif numMscans == 4 % VISTA 2 FBFB FBFB FB... = 0101 0101 23...
%        
%         for FrameNum = 2:2:size(vol,3)
%                 vol_reorder(:,:,FrameNum) = fliplr(vol(:,:,FrameNum)); % flip backward scans
%         end
%         
%         j=1;
%         k=4;
%         for i=1:numMscans*2:size(vol,3)-4  
%             vol_reorder(:,:,i:i+3) = vol(:,:,j:j+3);
%             j=j+4;
%         end
%         for i=numMscans+1:numMscans*2:size(vol,3)-4
%             vol_reorder(:,:,i:i+3) = vol(:,:,k+1:k+3);
%             k=k+4;
%         end
%     end
end


