function genDOPU_aviMontage(kernel,adapt,savepath,fn_num,...
    avgOCT,DOPU_test,DOPU_Bscans,cmap_dopu_r)

    % create movie file %
    % create the video writer with 1 fps
    knum=sprintf('%d%02d',kernel(1),kernel(2));
    if adapt==1
        vidname=fullfile(savepath,[fn_num,'adaptive',knum,'_montage.avi']);
    else
        vidname=fullfile(savepath,[fn_num,'rigid',knum,'_montage.avi']);
    end
    writerObj = VideoWriter(vidname);
    writerObj.FrameRate = 10; % set the seconds per image

    % open the video writer
    open(writerObj);
    for ii = 1:size(avgOCT,3)-1
        figure(1)  

        % logscale avgOCT
         if ii==1 || ii==size(avgOCT,3)
            G = flipud(squeeze(avgOCT(:,:,ii)));
        else
            % average one frame prior and after for display
            G = flipud(mean(avgOCT(:,:,ii-1:ii+1),3));
         end
        G(G==0)=30;
        imgL=10*log10(G); 
        imgL=mat2gray(imgL);
        imgL_test=imadjust(imgL,[0.25,0.9],[]);

        % OCTA
    %     imgA = mat2gray(flipud(OCTA(:,:,ii)));
    %     imgA = imadjust(imgA,[0,0.2],[],2.2);
        imgA=[];

        % DOPU
        imgD=flipud(DOPU_test(:,:,ii));

        % DOPU B-scans
        imgD_rgb_r = DOPU_Bscans(:,:,:,ii);
        [imgD_ind, cmap]=rgb2ind(imgD_rgb_r,32); % 64 is the highest number allowable if the bg is black

        ax1=subplot(3,1,1);imagesc(imgL_test);colormap(ax1,gray);title('log avgOCT');xlabel(ii);colorbar
    %     ax2=subplot(1,3,2);imagesc(imgA);colormap(ax2,gray);title('OCTA');xlabel(ii);
        ax2=subplot(3,1,2);imagesc(imgD_ind);colormap(ax2,cmap);title('DOPU B-scan');xlabel(ii);colorbar
        ax3=subplot(3,1,3);imagesc(imgD,[0,1]);colormap(ax3,cmap_dopu_r);colorbar;title('DOPU');xlabel(ii);
        set(gcf,'Position',[200 100 400 700])

        % convert the image to a frame
        frame = getframe(gcf) ;

        % write the frames to the video
        writeVideo(writerObj, frame);
        pause(0.01)
    end
    close(writerObj);

end