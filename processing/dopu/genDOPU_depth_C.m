function [DOPU_depth_C,DOPU_depth_C_bot]= genDOPU_depth_C(DOPU_test)
% genDOPU_depth_C estimates low melanin region boundaries based on DOPU
% information
%
% Input: Filtered DOPU volume DOPU_test (NEEDS to be both filtered and
% right-side up for this process!) 
% 
% % % example code for preparing DOPU
% % disp('Filtering DOPU...');
% % DOPU_filt(DOPU>0.95) = 1; % threshold theDOPU
% % DOPU_test=(medfilt3(DOPU_filt, [3 5 3])); % [depth width frames]
% % disp('DOPU filtered.');
% % 
% % DOPU_test=flipud(DOPU_test); % Need it right-side up
% 
% Outputs:  DOPU_depth_C - depth profile of RPE surface across volume
%           DOPU_depth_C_bot - depth profile of lower choroid boundary
%
% Note: some parameter adjustment may be necessary for improved quality.
% Data display/check code included (commented out) at bottom of function.

% size setup
numPoints=size(DOPU_test,1);
numAlines=size(DOPU_test,2);
numBscans=size(DOPU_test,3);

% prepare the martices for the segmentation of the DOPU and set to zeros
DOPU_depth_C = zeros(numAlines,numBscans); % RPE surface
DOPU_depth_C_bot = zeros(numAlines,numBscans); % choroid bottom boundary


%% Finding the cut lines (segmentation lines) of the low DOPU region
startZ = 1; %depth to start searching at (usd to cut out noise above retina)
thresh = 0.95; %threshold for finding low DOPU values (typically 0.95, lower values get closer to teh low dopur region)
smoothing = 10;% if noisy, do 30, else do 5-20
botcrop = 30; % crop out any FPN near bottom

for B=1:numBscans % loop across the frames
    for A=1:numAlines-1 % loop across the A-lines
        
        % extract A-line (all depth values)
        Aline=squeeze(DOPU_test(:,A,B));
        
        % smooth out any outliers
        Aline=smoothdata(Aline,'gaussian',smoothing); 
        
        % find all indices where DOPU < thresh
        idx = find(Aline(startZ:end-botcrop)<thresh)+startZ;
        
        % extract the top and bottom values (if any are found), else set
        % the DOPU depth at the bottom of the volume
        if idx
            DOPU_depth_C(A,B)=idx(1);
            DOPU_depth_C_bot(A,B)=idx(end);
        else
            DOPU_depth_C(A,B)=numPoints;
            DOPU_depth_C_bot(A,B)=numPoints;
        end
    end
    
    % copy the second-last A-line just to avoid zero values at edge
    DOPU_depth_C(end,:)=DOPU_depth_C(end-1,:);
    DOPU_depth_C_bot(end,:)=DOPU_depth_C_bot(end-1,:);
    
    % smooth out the final segmentation lines
    DOPU_depth_C(:,B)=smoothdata(DOPU_depth_C(:,B),'gaussian',50);
    DOPU_depth_C_bot(:,B)=smoothdata(DOPU_depth_C_bot(:,B),'gaussian',50);
    
%     fprintf('DOPU frame: %d\n', B);
end




%% check DOPU image
% % displaying an elevation map
% imgDDC = notchfilter(imcenhance(mat2gray((DOPU_depth_C))));%imcenhance (can do log10 on the dopu depth?
% 
% % display with inverted values(smaller numbers are higher elevation originally)
% figure;i1=imagesc(1-imgDDC,[0,1]);title('DOPU_depth, method C','Interpreter', 'none');
% colormap(hot);colorbar; axis equal; axis tight;
% 
% figure('pos',[50 50 1500 500]);
% for ii=1:10:size(DOPU_test,3)-1 %cplxOCTA, DOPU_test
% % for ii=500
%     ax1=subplot(1,2,1);
%     imgA = squeeze(DOPU_test(:,:,ii));
%     imagesc(imgA,[0,1]);colormap(ax1,cmap_dopu_r);colorbar;%colormap(gray, cmap_dopu);
%     title(ii);
%     hold on;
%     plot(DOPU_depth_C(:,ii),'black','LineWidth', 2);
%     plot(DOPU_depth_C_bot(:,ii),'black','LineWidth', 2);
%     hold off;
% 
%     pause(0.0001);
% %     colorbar;
% end
end