function Montage_movie_func_v2(savepath,filename,avgOCT,DOPU_Bscan,entropy)
    %% Load & Format Data
%     folder = 'G:\EDUADO\RP005_2022.05.01_OU\OD'; %directory where OCT volume is
%     cd(savepath);
%     addpath('C:\Users\coil_\Desktop\Jordan\jordan');
    close all
%     fname = '10_09_54-';
%     fname_dopu = [fname '_DOPU_Bscans.tiff'];
    %fname_dopu = [fname ,'-cplxOCT__DOPU_Bscans.tiff'];
%     fname_avgOCT = [fname '_avgOCT']; %file name of OCT volume
    %fname_avgOCT = [fname, '-cplxOCT__avgOCT'];

%     Data = load([fname_avgOCT,'.mat']); %load OCT volume
%     avgOCT = Data.avgOCT; %pull data out of structure
    avgOCT_grey = mat2gray(avgOCT); %grey scale image
    imsize = 500;
    skip = 1;
    previewColumn =8;
    
%     info = imfinfo(fname_dopu);
%     numframes = length(info);
    numframes = size(avgOCT,3);
    
%     for j = 1:numframes
%         frames(:,:,:,j) = imread(fname_dopu, j);
%     end
    
    % Enface
    enface = flip(squeeze(sum(avgOCT,1)),2);
    enface2 = (enface-min(min(enface)))/(max(max(enface))-min(min(enface))) * 256;
    enface2 = flip(flip(imresize(mat2gray(enface2),[imsize,imsize]),1),2);
    
    % Create Movie
    
    movie_name = fullfile(savepath,[filename,'_mov', '.avi']);
    writerObj = VideoWriter(movie_name);
    open(writerObj);
    
    n = 1;
    for i = 2:skip:numframes-2
    
        if i==1
            I = avgOCT_grey(:,:,1);
        elseif i==numframes
            I = avgOCT_grey(:,:,numframes);
        else
            I = mean(avgOCT_grey(:,:,i),3);
        end
    
        OCT_image = imadjust(mat2gray(20*log10(I)));
        OCT_image = flip(OCT_image,1);
        OCT_image = imresize(OCT_image, [imsize,imsize]);
        entropy_im = imresize(entropy(:,:,i), [imsize,imsize]);
    
        current_frame = DOPU_Bscan(:,:,:,i); %current frame from tiff stack
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
        montage({ enface3, OCT_image, current_frame,entropy_im}, 'Size', [1 4],'ThumbnailSize',[]);
        F = getframe(gcf);
        writeVideo(writerObj, F);

        % create preview montage
        if mod(i,round(numframes/(previewColumn+1)))==0
            ImMontage{n} = enface3;
            ImMontage{n+1} = OCT_image;
            ImMontage{n+2} = current_frame;
            ImMontage{n+3} = entropy_im;
            n = n+4;
        end
        
    
    end
    close(writerObj);
    
    for j = 0:3
        ImMontageReorder{1+j} = ImMontage{1+j};
        ImMontageReorder{5+j} = ImMontage{17+j};
        ImMontageReorder{9+j} = ImMontage{5+j};
        ImMontageReorder{13+j} = ImMontage{21+j};
        ImMontageReorder{17+j} = ImMontage{9+j};
        ImMontageReorder{21+j} = ImMontage{25+j};
        ImMontageReorder{25+j} = ImMontage{13+j};
        ImMontageReorder{29+j} = ImMontage{29+j};
    end

    close all
    figure(1),montage(ImMontageReorder, 'Size', [4 8],'ThumbnailSize',[]);
    saveas(gcf,fullfile(savepath,[patient_id,'_preview_dopu.tif']));
    close all

end