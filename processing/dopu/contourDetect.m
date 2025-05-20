function contour = contourDetect(input)
    % input = 2D image to have the contour detected. Should be right-side
    % up (contour to be detected shoudl be on the top of the image)
    %
    % contour = 1D array of depth values

%     input=flipud(input);

    % set up sizes
    numPoints=size(input,1);
    numAlines=size(input,2);

    % binarize the input
    input=imbinarize(mat2gray(input));
    % figure;imagesc(input);

    startZ = 1; %depth to start searching at (usd to cut out noise)
    thresh = 0.5; %threshold for finding low DOPU values
    botcrop = 0; % crop out any FPN near bottom, increase if there is a lot

    % prepare the martices for the segmentation of the DOPU and set to zeros)
    temp_contour = zeros(1,numAlines);

    % segment the RPE
    for A=1:numAlines-1 % loop across the A-lines

        % extract A-line (all depth values)
        Aline=squeeze(input(:,A));

        % find all indices where DOPU > thresh
        idx = find(Aline(startZ:end-botcrop)>thresh)+startZ;

        % extract the top  values (if any are found), else set
        % the DOPU depth at the bottom of the volume
        if idx
            temp_contour(A)=idx(1);
        else
            temp_contour(A)=numPoints;
        end

        % copy the second-last A-line just to avoid zero values
        temp_contour(end)=temp_contour(end-1);
    end



%     % check DOPU image
%     botlim=5;
%     lowlim =10;
%     uplim = 20;
% 
%     figure('pos',[50 50 1200 500]);
%     imgA = input;
%     imagesc(imgA,[0,1]);colorbar;%colormap(gray, cmap_dopu);
%     hold on;
%     plot(temp_contour,'red','LineWidth', 1);
%     hold off;
    
    % remove any areas that are clipping
    depth=zeros(size(temp_contour));
    idx = find(temp_contour < numPoints);
    depth=temp_contour(idx);
    % figure;plot(idx, depth, '.-');

    % find and remove outliers
    [outliers,out_x] = rmoutliers(depth,'movmedian',5,'SamplePoints',idx);
    % figure;plot(idx,depth,'b.-',idx(~TF),B,'r-')
    % legend('Input Data','Output Data')

    % interpolate/extrapolate any remaining points
    x_points=1:numAlines;
    interpolated = interp1(idx(~out_x),outliers,x_points,'linear','extrap');
%     figure;
%     plot(idx(~out_x),outliers,'.',x_points,interpolated,'-');
%     % xlim([0 2*pi]);
%     title('linear Interpolation');

    % smooth the data slightly
    interp_s=smoothdata(interpolated,'gaussian',15); 
%     figure;plot(x_points,interpolated,'-');
%     hold on;
%     plot(x_points,interp_s,'-');
%     hold off;


    % remove any clipping
    contour = interp_s;
    contour(contour<0) = 0;
    
%     % check final contour
%     figure;imagesc(input);
%     hold on;
%     plot(x_points,contour,'r-');
%     hold off;


end