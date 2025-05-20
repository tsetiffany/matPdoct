function save_pdoct_output2(OCT3d, DOP3d, file_dir,fullfile_name, cmap_dopu_r)
% Preset parameters
vert_scan       = 0;
bool_save_nii   = 1;

OCT3d = flip(flip(OCT3d,1),2);
DOP3d = flip(flip(DOP3d,1),2);
%%
% Load raw tiff 
% cd(file_dir)
% try
    [~,file_name,ext] = fileparts(fullfile_name);
    
%     switch ext
%         case'.tif'
%             disp('Loading tiff files..')
%             file_name = file_name(1:end-5);
%             if vert_scan==1
%                 OCT3d = FastTiff(fullfile(file_dir,[file_name,'_octv.tif']));
%                 DOP3d = FastTiff(fullfile(file_dir,[file_name,'_dopu.tif']))/(2^16-1);
%             else
%                 OCT3d = flip(FastTiff(fullfile(file_dir,[file_name,'_octv.tif'])),3);
%                 DOP3d = flip(FastTiff(fullfile(file_dir,[file_name,'_dopu.tif']))/(2^16-1),3);
%             end
%             file_name = file_name(end-10:end-5);
%         case '.mat'
%             disp('Loading mat files..')
%             file_name = file_name(1:end-4);
%             if vert_scan==1
%                 OCT3d = flip(load(fullfile(file_dir,[file_name,'_oct.mat'])).avgOCT,1);
%                 DOP3d = flip(load(fullfile(file_dir,[file_name,'_dopu.mat'])).DOPU,1);
%             else
%                 OCT3d = flip(flip(load(fullfile(file_dir,[file_name,'_oct.mat'])).avgOCT,1),3);
%                 DOP3d = flip(flip(load(fullfile(file_dir,[file_name,'_dopu.mat'])).DOPU,1),3);
%             end
%         otherwise
%             warning('Error: invalid file extension.')
%     end
    
    %% preview OCT
    % lin_Data = OCT3d;
    % log_data = 20*log10(OCT3d);
    % log_data(log_data<0)=0;
    % figure(1),
    % for i=1:10:size(OCT3d,3)
    %     imagesc(imadjust(mat2gray(log_data(:,:,i)))),colormap('gray')
    %     pause(0.1),
    % end
    % imagesc(abs(fft2(log_data(:,:,1)))),colormap('gray')
    
    %% DOPU filtering %%
    DOPU_filt=DOP3d;
    DOPU_filt(DOPU_filt>0.95) = 1; % threshold the DOPU
%     DOPU_filt=medfilt3(DOPU_filt,[3 5 3]); %[3 5(or 3) 3] (height, wdith, depth) with B-scans, [3 5 13] without B-scans
    
    %% Motion correction %%
    thresh = 50; % orig 5
    usfac = 1;
    numFrames = size(OCT3d, 3);
    numBatch = 2;
    %%% Set the shifting variable to save and analyze %%%
    xShift = zeros([numFrames 1]);
    yShift = zeros([numFrames 1]);
    
    OCT_mcorr = OCT3d;
    DOPU_mcorr = DOPU_filt;
    
    % Lateral registration
    %%% Start from third frame because first frame is not good for setting as reference frame %%%
%     ref_frame = round(numBatch / 2);
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
    
%     figure(2),
%     p = plot(yShift), title('y-shift'),
%     ylim([-thresh thresh])
%     exportgraphics(p,'motion_correction.png')
%     % pause(1),
%     close all
    % subplot(1,2,2),plot(yShift), title('y-shift')
    
    %% Frame averaging 
    frame_average = 2;
    
    linear_avgOCT = zeros(size(OCT3d,1),size(OCT3d,2),round(numFrames/numBatch));
    log_avgOCT = zeros(size(OCT3d,1),size(OCT3d,2),round(numFrames/numBatch));
    avgDOPU = zeros(size(OCT3d,1),size(OCT3d,2),round(numFrames/numBatch));
    
    for i = 1:numBatch:numFrames
        idx = ceil(i/numBatch);
        linear_avgOCT(:,:,idx) = mean(OCT_mcorr(:,:,i:i+frame_average-1),3);
        log_avgOCT(:,:,idx) = 20*log10(mean(OCT_mcorr(:,:,i:i+frame_average-1),3));
        avgDOPU(:,:,idx) = mean(DOPU_mcorr(:,:,i:i+frame_average-1),3);
    end
    log_avgOCT(log_avgOCT<0)=0;
    
    % generate composite DOPU
    DOPU_Bscans = genDOPU_combinedBscans(avgDOPU,linear_avgOCT,cmap_dopu_r);
    
    %% Create En face OCT stack
    
    enface_OCT_stack = zeros(size(OCT3d,2),size(OCT3d,3),3,round(numFrames/numBatch));
    enface_OCT =mat2gray(squeeze(mean(OCT3d,1)));
    enface_OCT_idx = gray2ind(enface_OCT,256);
    enface_OCT_rgb=ind2rgb(enface_OCT_idx,gray(256));
    
    for i = 1:numBatch:numFrames
        idx = ceil(i/numBatch);
        img_RGB = enface_OCT_rgb;
        if i == 1
            img_RGB(:,i:i+7,1)=255;
        elseif i == numFrames
            img_RGB(:,i-7:i,1)=255;
        else
            img_RGB(:,i-3:i+3,1)=255;
        end
    %     imInd = gray2ind(imgF,256);
    %     img_RGB=ind2rgb(imInd,hot(256));
        enface_OCT_stack(:,:,:,idx) = fliplr(img_RGB);
    end
    
    if vert_scan == 1
        exportTiff(enface_OCT_stack, fullfile(file_dir,[file_name,'_enface_oct_postion']))
    else
        exportTiff(imrotate(enface_OCT_stack,-90), fullfile(file_dir,[file_name,'_enface_oct_postion']))
    end
    
    %% Save
    % Save OCT B-scans %
    disp("- saving 3d oct ...")
    % save(fullfile(savepath,[patient_id,'_raw_oct.mat']), 'avgOCT', '-v7.3');
    exportTiff(fliplr(linear_avgOCT), fullfile(file_dir,[file_name,'_linear_oct']),1)
    exportTiff(fliplr(log_avgOCT), fullfile(file_dir,[file_name,'_log_oct']),1)
    % exportTiff(linear_avgOCT, fullfile(file_dir,[file_name,'_linear_oct']),1)
    % exportTiff(log_avgOCT, fullfile(file_dir,[file_name,'_log_oct']),1)
    
    % Save En face OCT image% 
    % disp("- saving en face ...")
    % f=flipud(imrotate(squeeze(mean(OCT3d,1)),-90));
    % % f=fliplr(squeeze(mean(OCT3d,1)));
    % f(f==0)=30;
    % imgF=10*log10(f); 
    % imgF=mat2gray(imgF);
    % imgF_test=imadjust(imgF,[0.25,0.9],[]);
    % figure(3);imshow((imgF));colormap(gray);axis equal; axis tight;
    % imwrite(imgF,fullfile(file_dir,[file_name,'_enface_oct.tif']));
    % close all
    
    % Save DOPU B-scans %
    disp("- saving 3d dopu ...")
    if vert_scan == 1
        exportTiff(flipud(DOPU_Bscans), fullfile(file_dir,[file_name,'_dopu_bscan']))
    else
        exportTiff(rot90(DOPU_Bscans,2), fullfile(file_dir,[file_name,'_dopu_bscan']))
    end
    
    % Save En face OCT image % 
    disp("- saving en face ...")
    if vert_scan == 1
        f=fliplr(squeeze(mean(OCT3d,1)));
    else
        f=flipud(imrotate(squeeze(mean(OCT3d,1)),-90));
    end
    
    f(f==0)=30;
    imgF=10*log10(f); 
    imgF=mat2gray(imgF);
%     imgF_test=imadjust(imgF,[0.25,0.9],[]);
%     figure(3);imshow((imgF));colormap(gray);axis equal; axis tight;
    imwrite(imgF,fullfile(file_dir,[file_name,'_enface_oct.tif']));
    close all
    
    % Save En face DOPU image % 
    disp("- saving en face ...")
    f=1-(squeeze(min(DOPU_filt,[],1)));
    % f(f==0)=30;
    if vert_scan == 1
        imgF=fliplr(f); 
    else
        imgF=imrotate(fliplr(f),-90); 
    end
    imgF=mat2gray(imgF);
    imInd = gray2ind(imgF,256);
    img_RGB=ind2rgb(imInd,hot(256));
    % imgF_test=imadjust(imgF,[0.25,0.9],[]);
%     figure(3);imshow((imgF));colormap(flipud(hot));axis equal; axis tight;
    % pause(0.1),
    % cmap = hot(64);
    imwrite(img_RGB,fullfile(file_dir,[file_name,'_enface_dopu.tif']));
    close all

        %% Create En face DOPU stack
    
    enface_stack = zeros(size(DOPU_filt,2),size(DOPU_filt,3),3,round(numFrames/numBatch));
    enface_im = 1-(squeeze(min(DOPU_filt,[],1)));
    enface_im_idx = gray2ind(enface_im,256);
    enface_im_rgb=ind2rgb(enface_im_idx,hot(256));
    
    for i = 1:numBatch:numFrames
        idx = ceil(i/numBatch);
        img_RGB = enface_im_rgb;
        if i == 1
            img_RGB(:,i:i+7,1)=255;
        elseif i == numFrames
            img_RGB(:,i-7:i,1)=255;
        else
            img_RGB(:,i-3:i+3,1)=255;
        end
    %     imInd = gray2ind(imgF,256);
    %     img_RGB=ind2rgb(imInd,hot(256));
        enface_stack(:,:,:,idx) = fliplr(img_RGB);
    end
    
    if vert_scan == 1
        exportTiff(enface_stack, fullfile(file_dir,[file_name,'_enface_dopu_postion']))
    else
        exportTiff(imrotate(enface_stack,-90), fullfile(file_dir,[file_name,'_enface_dopu_postion']))
    end
    
    % Save nifti %
    if bool_save_nii
    
        disp("- saving nifiti ...")
    
        % averaged OCT 
        raw = linear_avgOCT;
        raw_norm = (raw-min(raw(:)))/(max(raw(:)-min(raw(:))));
        raw_imadjust = (imadjustn(raw_norm));
        oct_uint8 = uint8(255*raw_imadjust);
        niftiwrite(permute(oct_uint8,[2,1,3]),fullfile(file_dir,[file_name,'_avgOCT']),'Compressed',true);
    
        % full OCT
%         dopu_uint8 = uint8(255*DOP3d);
%         niftiwrite(permute(dopu_uint8,[2,1,3]),fullfile(file_dir,[file_name(1:19),'_dopu']),'Compressed',false);
        raw = OCT3d;
        raw_norm = (raw-min(raw(:)))/(max(raw(:)-min(raw(:))));
        raw_imadjust = (imadjustn(raw_norm));
        oct_uint8 = uint8(255*raw_imadjust);
        niftiwrite(permute(oct_uint8,[2,1,3]),fullfile(file_dir,[file_name,'_OCT']),'Compressed',true);
    
    end

% catch
%     disp(['Error: file ' file_name 'not found.'])
% end