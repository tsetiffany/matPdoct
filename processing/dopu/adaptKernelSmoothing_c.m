function [out, adapt_angles, adapt_sizes]   =  adaptKernelSmoothing_c(input, kernel, depth)
    % input = input image (2D x 4 stokes vectors)
    % kernel = [z,x] averaging kernel, first number is axial, 2nd lateral
    % depth = 1D array of depths for the 2D image (forms a line)
    % out = output images after adaptive kernel smoothing (2D x 4 stokes
    % vectors)
    %
    % NOTE: input and depth must be same orientation!
%     figure;imagesc(input(:,:,1));hold on;plot(depth);hold off;

    imgDepth = size(input,1); % numPoints
    imgWidth = size(input,2); % numAlines
    
    % Get the angles
    h=(1:imgWidth)';
    dx = gradient(h, h);                               % Numerical Derivative: dx
    dy = gradient(depth, h);                               % Numerical Derivative: dy
    angle = -atan2d(dy, dx);                            % Angle Calculation (Degrees) 
    %(take negative because image display is vertically flipped from plot
    %display)

%     figure;imagesc(input);
%     hold on;plot(h,depth,'-');
%     angstr = sprintfc('\\angle %.0f \\circ',angle);
%     text(h, depth, angstr, 'HorizontalAlignment','left')
%     hold off;
%     
    

    out = zeros(size(input)); %prepare the output matrix the same size as the input frame

    %% loop setup and execution
    adapt_sizes = zeros(imgWidth,1);

    kern=ones([kernel(1),kernel(2)]); % set up kernel 
    
    for ii=1:imgWidth % loop across the scan
        rotkern = imrotate(kern, angle(ii)); % rotated kernel for this A-line angle
        rotkern=repmat(rotkern,[1,1,4]); % duplicate this for every Stokes vector
        
        % identify middle elements along width and depth
        kernW=size(rotkern,2);
        Mw=round(kernW./2);
        kernH=size(rotkern,1);
        Mh=round(kernH./2);
        
        % determine bounds for width of input slice
        cL = ii-(Mw-1); % left bound, subtract from current element
        if cL < 1 % check for border conditions
            %shift kernel start up number of positions overlapped
            rotkern=rotkern(:,(1+(abs(cL)+1)):end,:);
            cL=1; 
        end
        cR = ii+(kernW-Mw); % right bound, add from current element
        if cR > imgWidth
            cR= imgWidth;
            % set kernel end to width of slice
            rotkern=rotkern(:,1:(cR-cL+1),:);
        end 
        sliceW=input(:,cL:cR,:); % Extract width slice (speeds up execution)

        % Loop for first few rows
        for jj=1:Mh-1 
            cropkern=rotkern((1+(abs((jj-(Mh-1)))+1)):end,:,:);
            cT=1;
            cB = jj+(kernH-Mh); % bottom bound, add from current element

            temp=sliceW(cT:cB,:,:).*cropkern; % multiply the slices element by element
            
            out(jj,ii,:)=sum(temp,[1,2]); % sum all values for each vector
        end
        
        % Loop for middle depths    
        for jj=Mh:imgDepth-(kernH-Mh+1)        
            % determine bounds for width of input slice
            cT = jj-(Mh-1); % top bound, subtract from current element
            cB = jj+(kernH-Mh); % bottom bound, add from current element
            temp=sliceW(cT:cB,:,:).*rotkern; % multiply the slices element by element
            out(jj,ii,:)=sum(temp,[1,2]); % sum all values       
        end
        
        % Loop for last few rows
        for jj=imgDepth-(kernH-Mh+1)+1:imgDepth 
            cT = jj-(Mh-1); % top bound, subtract from current element
            cB = imgDepth;
            cropkern=rotkern(1:(cB-cT+1),:,:);
            
            temp=sliceW(cT:cB,:,:).*cropkern; % multiply the slices element by element
            out(jj,ii,:)=sum(temp,[1,2]); % sum all values
        end
        
        
        
        %%%%%%%%%%% Old version of code with one loop. Runs a tiny bit
        %%%%%%%%%%% slower so is commented out and left for posterity
%         for jj=1:imgDepth % loop down the depth
%             cropkern=rotkern; % set up cropped kernel for boundary conditions
%             
%             % determine bounds for width of input slice
%             cT = jj-(Mh-1); % top bound, subtract from current element
%             if cT <1 % check for border conditions
%                 %shift kernel start up number of positions overlapped
%                 cropkern=cropkern((1+(abs(cT)+1)):end,:);
%                 cT=1;
%             end
%             cB = jj+(size(rotkern,1)-Mh); % bottom bound, add from current element
%             if cB > imgDepth
%                 cB = imgDepth;
%                 % set kernel end to height of slice
%                 cropkern=cropkern(1:(cB-cT+1),:);
%             end
% 
%             temp=sliceW(cT:cB,:).*cropkern; % multiply the slices element by element
%             out(jj,ii)=sum(temp(:)); % sum all values
%             
%         end


        adapt_sizes(ii)=sum(rotkern(:,:,1),[1,2])-sum(kern(:)); % save the size difference from expected
    end


    % save angle and size variables
    adapt_angles = angle;
end