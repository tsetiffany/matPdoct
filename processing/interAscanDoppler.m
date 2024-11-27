function dVolume = interAscanDoppler(Volume, numMscans)

dVolume = zeros(size(Volume,1),size(Volume,2)/2,size(Volume,3)/numMscans);
dBscans = zeros(size(Volume,1),size(Volume,2)/2,numMscans);

for ima = 1:numMscans:size(Volume,3)
    K = ceil(ima/numMscans);
    
    for n = 0:numMscans-1
        
        Bscan = Volume(:,:,ima+n);
        
%         % A-line bulk phase correction %
%         for j=1:size(Bscan,2)-1
%             phi = angle(sum(Bscan(:,j+1).*conj(Bscan(:,j))));
%             Bscan(:,j+1) = Bscan(:,j+1) .* exp(-1j*phi);
%         end
        
        % inter Aline Doppler %
        for i=1:2:size(Bscan,2)
            dOCT(:,ceil(i/2)) = angle(Bscan(:,i).*conj(Bscan(:,i+1)));
        end
        
        % differentiate scan direction
        if mod(ima,numMscans*2) == 0 || mod(ima,numMscans*2) > numMscans
            dBscans(:,:,n+1) = -dOCT;
        else
            dBscans(:,:,n+1) = dOCT;
        end
        
    end
    dVolume(:,:,K) = squeeze(mean(dBscans,3));
end
    