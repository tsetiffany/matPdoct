function out = contourCorrect(input, depth, rem)
% Input = input volume (sould be right-side up, same orientation as the
% depth map in depth
% depth = depth map of contour to correct
% rem = 1, remove contour, 0 add back contour
% NOTE: input and depth must be same orientation!
    % figure;imagesc(input(:,:,100));hold on;plot(depth(:,100));hold off;
    
% setup output volume to be shifted
out=input;

% set sizes
% numPoints = size(input,1);
numAlines = size(input,2);
numFrames = size(input,3);

% set center of depth map as reference and shift everything relative to it
ref = depth(round(numAlines./2),round(numFrames./2));

% prepare the shift amount variable
if rem == 1
    % correcting the curve, make the shifts negative
    depth_shift = -(depth-ref);
else
    % putting curve back, make the shifts positive
    depth_shift = depth-ref;
end

% apply the shifts
for I = 1:numFrames
    for J = 1:numAlines
        out(:,  J, I) = circshift(input(:, J, I), [round(depth_shift(J,I)), 0]);
    end
end

end