function out = adaptKernelSmoothing_r(Bscan, kern, numChunks, add_angle)
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

%% parameter and input setup
% % test input parameters
% Bscan = S0_NC;
% % Bscan = S1_NC;
% % Bscan = S2;
% % Bscan = S3;
% kern = [3,5];
% numChunks = 5;

% get the inputs
input= Bscan;
kernel = kern;
% kernel = ones(kern(1),kern(2)); %unrotated kernel

imgDepth = size(input,1); % numPoints
imgWidth = size(input,2); % numAlines

% check if optional parameter was used
if isempty(add_angle)
   add_angle=0; % set default value for additional rotation angle 
end

step = round(imgWidth/numChunks); % dividing the B-scan into X equal chunks

% prepare weights
weight_Width_S = linspace(1,kernel(2)-1,kernel(2)-1)';
weight_Width_M = ones([(step+20)-2*(kernel(2)-1), 1])*kernel(2);
weight_Width_E = flipud(weight_Width_S);
weight_Width   = cat(1,weight_Width_S,weight_Width_M,weight_Width_E);
% e.g., kernel [1,8], input 372 pts by 1000 alines, weight_width =
% [1,2,3...8....8....3,2,1] for 1000 values.  Convolution kernels are all
% ones adn then weighted by the weight_width afterwards

weight_Depth_S = linspace(1,kernel(1)-1,kernel(1)-1)';
weight_Depth_M = ones([imgDepth-2*(kernel(1)-1), 1])*kernel(1);
weight_Depth_E = flipud(weight_Depth_S);
weight_Depth   = cat(1,weight_Depth_S,weight_Depth_M,weight_Depth_E);

out = zeros(size(input)); %prepare the output matrix the same size as the input frame

%% loop setup and execution
count = 1; % count how many slices have been processed

for ii=1:step:imgWidth % loop over the X chunks of the scan
%     ii=1;
    if count == 1
        % take an extra 20 pixels after
        slice = input(:,ii:ii+(step-1)+20);
    elseif count < numChunks
        % take an extra total 20 pixels after and beore
        slice = input(:,ii-10:ii+(step-1)+10);
    else
        % last slice, have to take previous 20
        slice = input(:,ii-20:end);
    end

%     %display the image (for testing purposes)
%     figure;
%     subplot(2,1,1);
%     imagesc(slice);
%     title('Bscan slice');
%     subplot(2,1,2);
%     imagesc(imbinarize(slice));
%     title('Binarized');
%     axis on, axis normal, hold on;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % use regionprops on the binarized image to find the orientation and
    % major axis length/area
    stats = regionprops(imbinarize(input), 'Orientation','MajorAxisLength','MinorAxisLength','Area'); %orientation returned wrt x axis from -90 to 90, positive = ccw, negative = cw

    % the largest major axis length/area correcpond, hopefully, to the main
    % RPE layer (neeed the cat 1 to get all the values)
    theta = cat(1,stats.Orientation);
    mal = cat(1,stats.MajorAxisLength);
    mil = cat(1,stats.MinorAxisLength);
    area = cat(1,stats.Area); % may not need the area, but keep this in
    %     case tests are inconclusive, may need area instead of Max axis length
    %     in cases of ONH

    % Extract the angle by findign the one for the longest continous region
    % (more than likely coresponds to the RPE)
    angle = theta(area==max(area(:)));

    % Find the major and minior axis lengths corresponding to the largest
    % continuous area (again most liketly the RPE)
    max_alen = mal(area==max(area(:)));
    min_alen = mil(area==max(area(:)));
    
    % add in any additional rotation requested
    tot_angle = angle + add_angle;
   
    
%     %%%%%%%%%%%%% Elliptical kernel %%%%%%%%%%%%%%%%%%%%
%     l=(kernel(2)+1)./2; %Length
%     w=(kernel(1)+1)./2; %Width
%     [X, Y] = meshgrid(-kernel(2):kernel(2),-kernel(1):kernel(1)); %make a meshgrid: use the size of your image instead
%     ellipse = (X/l).^2+(Y/w).^2<=1; %Your Binary Mask which you multiply to your image, but make sure you change the size of your mesh-grid
%     rotkern = double(imrotate(ellipse,tot_angle));
    
    
    
    
%     %%%%%%%%%%%%%%%%%%%%% averaging method %%%%%%%%%%%%%%%%%%%
%     % apply rotation and average
%     kern = ones(kernel(1), kernel(2));
%     rotkern = imrotate(kern, tot_angle); 
%     img_con = conv2(slice, rotkern, 'same'); % if rotated kernel is symmetrical, no need to add an extra rot90x2
%     outputX = img_con./floor(sum(rotkern(:)));
%    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%% start the convolutions %%%%%%%%%%%%%%%%%%%%%%
    outputY = zeros(size(slice));

    % apply rotation
    kernelY  = ones([1 kernel(2)]);
    rotkernY = imrotate(kernelY, tot_angle); 
    img_conY = conv2(slice, rotkernY, 'same'); % if rotated kernel is symmetrical, no need to add an extra rot90

    for k = 1:imgDepth
        outputY(k,:) = img_conY(k,:)'./weight_Width;
    end

    outputX = outputY;

    kernelX  = ones([kernel(1),1]);
    % appy rotation
    rotkernX = imrotate(kernelX, tot_angle);
    img_conX = conv2(outputY, rotkernX, 'same');

    for kk = 1:size(slice,2)
        outputX(:,kk) = img_conX(:,kk)./weight_Depth;
    end
     

    if count == 1
        % put the results into the output matrix (withtou last 20)
        out(:,ii:(ii+step-1)) = outputX(:,1:end-20);
    elseif count < numChunks
        % put the results into the output matrix
        out(:,ii:(ii+step-1)) = outputX(:,1+10:end-10);
    else
        % put the results into the output matrix (withtou first 20)
        out(:,ii:(ii+step-1)) = outputX(:,1+20:end);
    end
   
    count  = count+1;
end

%%
% 
% 
% 
% 
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