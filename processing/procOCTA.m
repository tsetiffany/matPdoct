function [avgOCT, OCTA, xShift, yShift] = procOCTA(cplxOCT, numMscans, typeOCTA)
%%% typeOCTA = 0 Var; 1 Sub; 2 phaseVariance
    
%%% 2D motion correction %%%
S_Mcorr  = cplxOCT;
numFrames = size(cplxOCT,3);

usfac  = 1;
xShift = zeros([numFrames 1]);
yShift = zeros([numFrames 1]);

for I = 1:numMscans:numFrames-1
	for j=1:numMscans-1

		[output, ~] = dftregistration(fft2(20.*log10(abs(S_Mcorr(:, :, I)))),...
					fft2(20.*log10(abs(S_Mcorr(:, :, I+j)))), usfac);

		xShift(I+j) = round(output(4)); 
		yShift(I+j) = round(output(3));

		S_Mcorr(:, :, I+j)  = circshift(S_Mcorr(:, :, I+j),  [round(output(3)) round(output(4))]);
		
	end
end

%% OCT Average & OCTA Processing %%
volume = S_Mcorr;

if numMscans == 2

	for I = 1:numMscans:size(volume,3)-1
		K = ((I-1)/numMscans)+1;
		
		Xconj_1   = volume(:,:,I+1).*conj(volume(:,:,I));
		BulkOff_1 = repmat(angle(sum(Xconj_1,1)), [size(Xconj_1,1) 1]);
		
		Bscan_1  = volume(:,:,I);
		Bscan_2  = volume(:,:,I+1).*exp(-1j*BulkOff_1);
		
		% Average OCT *
		avgOCT(:,:,K) = (abs(Bscan_1) + abs(Bscan_2))./2;
		
		% Variance %
		Var(:,:,K) = abs(var(cat(3,Bscan_1,Bscan_2),0,3));
		
		% Substraction %
		Sub(:,:,K) = abs(Bscan_1 - Bscan_2); 
		
		% Phase Variance %
		avgOCT_conj = conj(Bscan_1 + Bscan_2);
		phiVar(:,:,K) = abs(angle(Bscan_1.*avgOCT_conj)) + abs(angle(Bscan_2.*avgOCT_conj));
        
		fprintf('OCTA volume process: %d\n', K);
    end
	
elseif numMscans == 3

	for I = 1:numMscans:size(volume,3)
		K = ((I-1)/numMscans)+1;
		
		Xconj_2   = volume(:,:,I+1).*conj(volume(:,:,I));
		Xconj_3   = volume(:,:,I+2).*conj(volume(:,:,I));
		BulkOff_2 = repmat(angle(sum(Xconj_2)), [size(Xconj_2,1) 1]);
		BulkOff_3 = repmat(angle(sum(Xconj_3)), [size(Xconj_3,1) 1]);
		
		Bscan_1  = volume(:,:,I);
		Bscan_2  = volume(:,:,I+1) .* exp(-1j*BulkOff_2);
		Bscan_3  = volume(:,:,I+2) .* exp(-1j*BulkOff_3);
		
		% Average OCT %
		avgOCT(:,:,K) = (abs(Bscan_1) + abs(Bscan_2) + abs(Bscan_3))/3;
		
		% Variance %
		Var(:,:,K) = abs(var(cat(3,(Bscan_1),(Bscan_2),(Bscan_3)),0,3));
		
		% Phase Variance %
		avgOCT_conj = conj(Bscan_1 + Bscan_2 + Bscan_3);
		phiVar(:,:,K) = abs(angle(Bscan_1.*avgOCT_conj))...
			+ abs(angle(Bscan_2.*avgOCT_conj))...
			+ abs(angle(Bscan_3.*avgOCT_conj));
			
		% Substraction %
		Sub(:,:,K) = (abs(Bscan_1 - Bscan_2)+abs(Bscan_2 - Bscan_3)); 
		
%		% Decorrelation %
%		Dec(:,:,K) = 1 - ((abs(Bscan_1).*abs(Bscan_2))...
%			./((abs(Bscan_1).^2 + abs(Bscan_2).^2)./2));
		
		fprintf('OCTA volume process: %d\n', K);
    end    
    
elseif numMscans == 4

	for I = 1:numMscans:size(volume,3)
		K = ((I-1)/numMscans)+1;
		
		Xconj_2   = volume(:,:,I+1).*conj(volume(:,:,I));
		Xconj_3   = volume(:,:,I+2).*conj(volume(:,:,I));
		Xconj_4   = volume(:,:,I+3).*conj(volume(:,:,I));
		BulkOff_2 = repmat(angle(sum(Xconj_2)), [size(Xconj_2,1) 1]);
		BulkOff_3 = repmat(angle(sum(Xconj_3)), [size(Xconj_3,1) 1]);
		BulkOff_4 = repmat(angle(sum(Xconj_4)), [size(Xconj_4,1) 1]);
		
		Bscan_1  = volume(:,:,I);
		Bscan_2  = volume(:,:,I+1) .* exp(-1j*BulkOff_2);
		Bscan_3  = volume(:,:,I+2) .* exp(-1j*BulkOff_3);
		Bscan_4  = volume(:,:,I+3) .* exp(-1j*BulkOff_4);
		
		% Average OCT %
		avgOCT(:,:,K) = (abs(Bscan_1) + abs(Bscan_2) + abs(Bscan_3) + abs(Bscan_4))/4;
		
		% Variance %
		Var(:,:,K) = abs(var(cat(3,(Bscan_1),(Bscan_2),(Bscan_3),(Bscan_4)),0,3));
		
		% Phase Variance %
		avgOCT_conj = conj(Bscan_1 + Bscan_2 + Bscan_3 + Bscan_4);
		phiVar(:,:,K) = abs(angle(Bscan_1.*avgOCT_conj))...
			+ abs(angle(Bscan_2.*avgOCT_conj))...
			+ abs(angle(Bscan_3.*avgOCT_conj))...
			+ abs(angle(Bscan_4.*avgOCT_conj));
			
		% Substraction %
		Sub(:,:,K) = (abs(Bscan_1 - Bscan_2)+abs(Bscan_2 - Bscan_3)+abs(Bscan_3 - Bscan_4)); 
		
%		% Decorrelation %
%		Dec(:,:,K) = 1 - ((abs(Bscan_1).*abs(Bscan_2))...
%			./((abs(Bscan_1).^2 + abs(Bscan_2).^2)./2));
		
		fprintf('OCTA volume process: %d\n', K);
    end
    
elseif numMscans == 1
    
    avgOCT  = cplxOCT;
    Var     = abs(cplxOCT);
    Sub     = abs(cplxOCT);
    phiVar  = abs(cplxOCT);
    
elseif numMscans == 10 %test purpose
    
    	for I = 1:numMscans:size(volume,3)
            K = ((I-1)/numMscans)+1;
            
            Bscans = volume(:,:,I);
            for ii = 1:numMscans-1
                Xconj = volume(:,:,I+ii).*conj(volume(:,:,I));
                BulkOff = repmat(angle(sum(Xconj)), [size(Xconj,1) 1]);
                Bscans(:,:,ii+1)  = volume(:,:,I+ii) .* exp(-1j*BulkOff);
            end
            
            % Average OCT %
            avgOCT(:,:,K) = mean(abs(Bscans(:,:,1:2)),3);

            % Variance %
            Var(:,:,K) = abs(var(Bscans(:,:,1:2),0,3));

%             % Phase Variance %
%             avgOCT_conj = conj(Bscan_1 + Bscan_2 + Bscan_3 + Bscan_4);
%             phiVar(:,:,K) = abs(angle(Bscan_1.*avgOCT_conj))...
%                 + abs(angle(Bscan_2.*avgOCT_conj))...
%                 + abs(angle(Bscan_3.*avgOCT_conj))...
%                 + abs(angle(Bscan_4.*avgOCT_conj));
% 
            % Substraction %
            Sub(:,:,K) = abs(Bscans(:,:,1) - Bscans(:,:,6)); 

    %		% Decorrelation %
    %		Dec(:,:,K) = 1 - ((abs(Bscan_1).*abs(Bscan_2))...
    %			./((abs(Bscan_1).^2 + abs(Bscan_2).^2)./2));
    
            if mod(K,2)==0
                avgOCT(:,:,K) = fliplr(avgOCT(:,:,K));
                Var(:,:,K) = fliplr(Var(:,:,K));
                Sub(:,:,K) = fliplr(Sub(:,:,K));
            end
            
            fprintf('OCTA volume process: %d\n', K);
        end
    
end

switch typeOCTA
	case 'var'
		OCTA = Var;
	case 'sub'
		OCTA = Sub;
	case 'phi'
		OCTA = phiVar;
end		

