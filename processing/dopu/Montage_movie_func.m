function Montage_movie_func(folder, fname, avgOCT, imsize)
    %% Load & Format Data
%     folder = 'G:\EDUADO\RP005_2022.05.01_OU\OD'; %directory where OCT volume is
    cd(folder);
    fname_dopu = [fname '_DOPU_Bscans.tiff'];
    %fname_dopu = [fname ,'-cplxOCT__DOPU_Bscans.tiff'];
    fname_avgOCT = [fname '_avgOCT']; %file name of OCT volume
    %fname_avgOCT = [fname, '-cplxOCT__avgOCT'];

%     Data = load([fname_avgOCT,'.mat']); %load OCT volume
%     avgOCT = Data.avgOCT; %pull data out of structure
    avgOCT_grey = mat2gray(avgOCT); %grey scale image
%     imsize = 250;
    skip = 2;
%     videoType = 'Archival';
    
    info = imfinfo(fname_dopu);
    numframes = length(info);
    
    for j = 1:numframes
        frames(:,:,:,j) = imread(fname_dopu, j);
    end
    
    % Enface
    enface = flip(squeeze(sum(avgOCT,1)),2);
    enface2 = (enface-min(min(enface)))/(max(max(enface))-min(min(enface))) * 256;
    enface2 = flip(flip(imresize(mat2gray(enface2),[imsize,imsize]),1),2);
    
    % Create Movie
    
    movie_name = [fname,'_preview', '.avi'];
    writerObj = VideoWriter(movie_name,'Motion JPEG AVI');
    open(writerObj);
    
    
    for i = 2:skip:numframes-2
    
        if i==1
            I = avgOCT_grey(:,:,1);
        elseif i==numframes
            I = avgOCT_grey(:,:,numframes);
        else
            I = mean(avgOCT_grey(:,:,i-1:i+1),3);
        end
    
        OCT_image = imadjust(mat2gray(20*log10(I)));
        OCT_image = flip(OCT_image,1);
        OCT_image = imresize(OCT_image, [imsize,imsize]);
    
        current_frame = frames(:,:,:,i); %current frame from tiff stack
        current_frame = imresize(current_frame, [imsize,imsize]);
    
        figure(1)
        imshow(enface2)
        hold on
        p1 = [i*imsize/numframes, imsize];
        p2 = [i*imsize/numframes, 0];
        plot([p1(1),p2(1)], [p1(2),p2(2)], 'Color', 'r', 'LineWidth',2);
        set(gca,'position',[0 0 1 1],'units','normalized')
        F1 = getframe(gcf);
        [enface3, Map] = frame2im(F1);
        %enface3 = RemoveWhiteSpace(enface3);
        enface3 = imresize(enface3, [imsize, imsize]);
        hold off
    
        figure(2)
        montage({ enface3, OCT_image, current_frame}, 'Size', [1 3],'ThumbnailSize',[]);
        F = getframe(gcf);
        writeVideo(writerObj, F);
    
    end
    
    close(writerObj);


end