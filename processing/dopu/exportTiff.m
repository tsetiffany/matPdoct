function exportTiff(vol, fname, norm)

if nargin < 3
    norm = 0;   % default setting
end

% vol = avgOCT_mcorr;
% fname = [fname_save, '_avgOCT'];

if size(vol,4) ~= 1
    %%% save image as RGB tif stack %%% 
    t=Tiff([fname,'.tiff'],'w');
    tagstruct.ImageLength = size(vol,1); % image height
    tagstruct.ImageWidth = size(vol,2); % image width
    tagstruct.Photometric = Tiff.Photometric.RGB;
    % tagstruct.Photometric = Tiff.Photometric.LinearRaw; % https://de.mathworks.com/help/matlab/ref/tiff.html
    tagstruct.BitsPerSample = 8;
    tagstruct.SamplesPerPixel = 3;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky; % groups rgb values into a single pixel instead of saving each channel separately for a tiff image
    tagstruct.Software = 'MATLAB';
    setTag(t,tagstruct)
    write(t,squeeze(im2uint8(((vol(:,:,:,1))))));
    for i = 2:size(vol,4)
    %     ImageLayer = vol(:,:,i);
        writeDirectory(t);
        setTag(t,tagstruct)
        write(t,squeeze(im2uint8(((vol(:,:,:,i)))))) %%%appends the next layer to the same file t
    end
    close(t) 
else
    %%% save image as grayscale tif stack %%% 
    t=Tiff([fname,'.tiff'],'w');
    tagstruct.ImageLength = size(vol,1); % image height
    tagstruct.ImageWidth = size(vol,2); % image width
%     tagstruct.Photometric = Tiff.Photometric.RGB;
    tagstruct.Photometric = Tiff.Photometric.LinearRaw; % https://de.mathworks.com/help/matlab/ref/tiff.html
    tagstruct.BitsPerSample = 8;
    tagstruct.SamplesPerPixel = 1;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky; % groups rgb values into a single pixel instead of saving each channel separately for a tiff image
    tagstruct.Software = 'MATLAB';
    setTag(t,tagstruct)
    if norm == 0
        write(t,squeeze(im2uint8(((vol(:,:,1))))));
    else
        write(t,squeeze(im2uint8(imadjust(mat2gray(vol(:,:,1))))));
    end
    for i = 2:size(vol,3)
    %     ImageLayer = vol(:,:,i);
        writeDirectory(t);
        setTag(t,tagstruct)
        if norm == 0
            write(t,squeeze(im2uint8(((vol(:,:,i)))))) %%%appends the next layer to the same file t
        else
            write(t,squeeze(im2uint8(imadjust(mat2gray(vol(:,:,i)))))) %%%appends the next layer to the same file t
        end
    end
    close(t) 
end