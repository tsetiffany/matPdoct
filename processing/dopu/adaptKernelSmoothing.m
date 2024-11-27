function out = adaptKernelSmoothing(Bscan, kern, numChunks, add_angle)
%{
Adaptive kernel for RPE smoothing; Inputs: the PS frame to be proecssed, a
default kernel size, number of slices (should be 5), offset angle for rotation (optional if want to apply a general
rotation to the kernel)

Take the frame from the noise-corrected DOPU that is to be smoothed

Divide B-scan into five to ten pieces.  Easier, less time-consuming, and more
daptable to each volume size.

Find the angle of the RPE by using the orientation and major axis length in
regionprops.  Max axis length is the rpe, typically, and orientation is teh
angle

Create a kernel based on inputs.  Pad with zeros to allow for rotation.
Rotate using imrotate

Convolve/(find term for convolution without flipping) for all elements in
that chunk (and perhaps one or two more on the edge to allow for
overflow; need to revise tht for teh last chunk) with the rotated kernel.
Save the results into the output matrix.

Repeat for each A-line chunk.

%}

%% ALTERNATE:
% use the same code as the smooth2Dfilter, only rotate the diretional
% kernel used!

% % test input parameters
% Bscan = S0_NC;
% % Bscan = S1_NC;
% % Bscan = S2;
% % Bscan = S3;
% kern = [3,5];
% numChunks = 5;

input= Bscan;
kernel = kern;

imgDepth = size(input,1);
imgWidth = size(input,2);

weight_Width_S = linspace(1,kernel(2)-1,kernel(2)-1)';
weight_Width_M = ones([imgWidth-2*(kernel(2)-1), 1])*kernel(2);
weight_Width_E = flipud(weight_Width_S);
weight_Width   = cat(1,weight_Width_S,weight_Width_M,weight_Width_E);
% e.g., kernel [1,8], input 372 pts by 1000 alines, weight_width =
% [1,2,3...8....8....3,2,1] for 1000 values.  Convolution kernels are all
% ones adn then weighted by the weight_width afterwards



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use regionprops on the binarized image to find the orientation and
% major axis length/area
stats = regionprops(imbinarize(input), 'Orientation','MajorAxisLength','MinorAxisLength','Area'); %orientation returned wrt x axis from -90 to 90, positive = ccw, negative = cw

% the largest major axis length/area correcpond, hopefully, to the main
% RPE layer (neeed the cat 1 to get all the values)
theta = cat(1,stats.Orientation);
% mal = cat(1,stats.MajorAxisLength);
% mil = cat(1,stats.MinorAxisLength);
area = cat(1,stats.Area); % may not need the area, but keep this in
%     case tests are inconclusive, may need area instead of Max axis length
%     in cases of ONH

% Extract the angle by findign the one for the longest continous region
% (more than likely coresponds to the RPE)
angle = theta(area==max(area(:)));

% Find the major and minior axis lengths corresponding to the largest
% continuous area (again most liketly the RPE)
% max_alen = mal(area==max(area(:)));
% min_alen = mil(area==max(area(:)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



outputY = zeros(size(input));


% appy rotation
kernelY  = ones([1 kernel(2)]);
rotkernY = imrotate(kernelY, angle);
img_conY = conv2(input, rotkernY, 'same');

for k = 1:imgDepth
    outputY(k,:) = img_conY(k,:)'./weight_Width;
end

weight_Depth_S = linspace(1,kernel(1)-1,kernel(1)-1)';
weight_Depth_M = ones([imgDepth-2*(kernel(1)-1), 1])*kernel(1);
weight_Depth_E = flipud(weight_Depth_S);
weight_Depth   = cat(1,weight_Depth_S,weight_Depth_M,weight_Depth_E);

outputX = outputY;

kernelX  = ones([kernel(1),1]);
% appy rotation
rotkernX = imrotate(kernelX, angle);
img_conX = conv2(outputY, rotkernX, 'same');

for kk = 1:imgWidth
    outputX(:,kk) = img_conX(:,kk)./weight_Depth;
end

out = outputX;



% %% parameter and input setup
% 
% %%%%%%%%%%%%%% testing input parameters %%%%%%%%%%%%%%
% Bscan = S0_NC;
% % Bscan = S1_NC;
% % Bscan = S2;
% % Bscan = S3;
% % kern = [3,5];
% % numChunks = 5;
% add_angle=0;
% 
% % prepare parameters from given inputs
% % % check if optional parameter was
% % if isempty(add_angle)
% %    add_angle=0; % set default value for additional rotation angle 
% % end
% 
% % numPoints = size(Bscan,1);
% numAlines = size(Bscan,2);
% kernel = ones(kern(1),kern(2)); %unrotated kernel (repalce with input params)
% step = round(numAlines/numChunks); % dividing the B-scan into X equal chunks
% 
% out = zeros(size(Bscan)); %prepare the output matrix the same size as the input frame
% 
% 
% %% loop setup and execution
% count = 1; % count how many slices have been processed
% 
% for ii=1:step:numAlines % loop over the X chunks of the scan
% %     ii=1;
%     if count < 5
%         % take one extra chunk for convolution to avoid cutoffs at borders
%         slice = Bscan(:,ii:ii+(step*2-1));
%     else
%         % last slice, have to take previous slice instead
%         slice = Bscan(:,ii-step:end);
%     end
%     
% %     % find the location of the RPE at the beginning and end of teh edge (or
% %     % perhaps somehow fit a line to the RPE roughly and get the angle from
% %     % that?)
% %     edges = edge(imbinarize(slice),'Canny');
% %     
% %     % the Hough transform is used to detect lines in an image, with the theta
% %     % built in to boot (which is theta + 90 degrees, clockwise wrt positive
% %     % x-axis)
% %     [H, T, R] = hough(imbinarize(slice)); % hough, theta, rho
% %     
% %     % display the image
% %     figure;
% %     subplot(3,1,1);
% %     imagesc(slice);
% %     title('Bscan slice');
% %     subplot(3,1,2);
% %     imagesc(edges);
% %     title('Edges');
% %     subplot(3,1,3);
% %     imshow(imadjust(rescale(H)),'XData',T,'YData',R,...
% %           'InitialMagnification','fit');
% %     title('Hough transform');
% %     xlabel('\theta'), ylabel('\rho');
% %     axis on, axis normal, hold on;
% %     colormap(gca,hot);
% % %
% %     % Need to figure out which theta has the most tendency or seomthing.....
%     
% 
% 
% %     %display the image (for testing purposes)
% %     figure;
% %     subplot(2,1,1);
% %     imagesc(slice);
% %     title('Bscan slice');
% %     subplot(2,1,2);
% %     imagesc(imbinarize(slice));
% %     title('Binarized');
% %     axis on, axis normal, hold on;
% 
%     % use regionprops on the binarized image to find the orientation and
%     % major axis length/area
%     stats = regionprops(imbinarize(slice), 'Orientation','MajorAxisLength','Area'); %orientation returned wrt x axis from -90 to 90, positive = ccw, negative = cw
%     
%     % the largest major axis length/area correcpond, hopefully, to the main
%     % RPE layer (neeed the cat 1 to get all the values)
%     theta = cat(1,stats.Orientation);
%     mal = cat(1,stats.MajorAxisLength);
% %     area = cat(1,stats.Area); % may not need the area, but keep this in
% %     case tests are inconclusive
%     
%     % Extract the angle by findign teh one for the longest continous region
%     % (more than likely coresponds to the RPE)
%     angle = theta(mal==max(mal(:)));
%     
%     % add in any additional rotation requested
%     tot_angle = angle + add_angle;
%     % apply rotation to the kernel by said angle
%     rotkern = imrotate(kernel, tot_angle); % keepingother settings deault; nearest neighbour interpolation and loose to allow fo rlarger boudnign box and preserving all teh kernel data
%     
%     % change the kernel to a smoothing/averaging kernel (or Gaussian?
%     % Double-check the implementation in smoothing 2D)
%     % (seems to do each direction independently, with more weight/division in the middle and less at ends?)
%     rot_smooth_kern = rotkern./sum(rotkern(:)); 
%     rot_smooth_kern = rot90(rot_smooth_kern,2); % flip the kernel (this is the equivalent of flipud(fliplr()))
%     % ^ not sure if the above is necessary
%     
%     % convolve with the chunk + 1 (determiend at top of loop)
%     smooth_slice = conv2(slice,rot_smooth_kern,'same'); % keep the output size the same by cropping off the excess
%     
%     % check if we're at the last slice yet
%     if count < 5
%         % put the results into the output matrix
%         out(:,ii:(ii+step-1)) = smooth_slice(:,1:step);
%     else
%         % put the results into the output matrix (second half)
%         out(:,ii:(ii+step-1)) = smooth_slice(:,step+1:end);
%     end
%     
%     count  = count+1;
% end
% 
% %%
% % 
% % 
% % 
% % 
% % normal smoothing with the 2D smohting filter
% nonAd = smooth2DFilter(Bscan,kern);
% 
% %display the images
% figure;
% subplot(4,1,1);
% imagesc(Bscan(1:120,:));title('Original bscan slice');
% subplot(4,1,2);
% imagesc(nonAd(1:120,:));title('Non adaptive smoothing kernel');
% subplot(4,1,3);
% imagesc(out(1:120,:));title('Adaptive averaging kernel');
% subplot(4,1,4);
% imagesc(out_smooth(1:120,:));title('Adaptive smoothing kernel');
% axis on, axis normal;
% colormap(hsv);



end