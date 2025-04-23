function [avgOCT, avgDOPU] = group_average_4(OCT3d, DOPU_filt)
    
    % Motion correction %%
    thresh = 50; % orig 5
    usfac = 1;
    numFrames = size(OCT3d, 3);
    numBatch = 4;
    %%% Set the shifting variable to save and analyze %%%
    xShift = zeros([numFrames 1]);
    yShift = zeros([numFrames 1]);
    frame_average = numBatch;
    
    
    % Lateral registration
    %%% Start from third frame because first frame is not good for setting as reference frame %%%
%     ref_frame = round(numBatch / 2); 
    if numBatch >1
        OCT_mcorr = OCT3d;
        DOPU_mcorr = DOPU_filt;
        for I = 1:numBatch:numFrames
	        for j= 1:numBatch-1
        
		        [output, ~] = dftregistration(fft2(OCT_mcorr(:, :, I)),...
					        fft2(OCT_mcorr(:, :, I+j)), usfac);
        
        % 		xShift(I+j) = round(output(4)); 
		        yShift(I+j) = round(output(3));
        
        
                %%% Thresholding  value was found via plotting the shifting values %%%
                if abs(output(3)) >= thresh
                    output(3) = 0;
                end
                
		        OCT_mcorr(:, :, I+j)  = circshift(OCT_mcorr(:, :, I+j),  [round(output(3)) 0]);
                DOPU_mcorr(:, :, I+j)  = circshift(DOPU_mcorr(:, :, I+j),  [round(output(3)) 0]);
        
            end
            if mod(I-1,100)==0
                fprintf('motion correction %d/%d \n',I-1,numFrames);
            end
        end

        % frame average
        for i = 1:numBatch:numFrames
            idx = ceil(i/numBatch);
            avgOCT(:,:,idx) = mean(OCT_mcorr(:,:,i:i+frame_average-1),3);
            avgDOPU(:,:,idx) = mean(DOPU_mcorr(:,:,i:i+frame_average-1),3);
        end

    else
        avgOCT = OCT3d;
        avgDOPU = DOPU_filt;
    end

end
