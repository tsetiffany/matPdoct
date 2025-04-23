function vol_bgsub = linear_intensity_bgsub(vol)

    enface_mask = zeros(size(vol,2),size(vol,3));
    enface_mask(2:31,2:31) = 1;
      
    
    noise_bg = zeros(size(vol,1),1);
    count = 0;
    for frame=1:size(vol,3)
        for aline=1:size(vol,2)
            if count >= 1000
                break
            end
            if enface_mask(aline,frame)==1
                noise_bg = noise_bg + vol(:,aline,frame);
                count = count + 1;
            end
        end
    end
    noise_bg = noise_bg / count;

    vol_bgsub = abs(vol - noise_bg);

%     figure(3),
%     imshow(log10(vol_bgsub(:,:,end-1)),[]);

%     figure(4),
%     subplot(1,3,1),plot(vol(:,:,end-1));
%     subplot(1,3,2),plot(vol_bgsub(:,:,end-1));
%     subplot(1,3,3),plot(noise_bg);

end