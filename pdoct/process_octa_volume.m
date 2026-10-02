function [avgOCT, subOCTA] = process_octa_volume(cplxData, numBMscans, usfac)
%PROCESS_OCTA_VOLUME Sub-pixel motion correction + bulk-phase-corrected
% OCTA calculation for a single channel's interleaved, cropped complex
% volume (4BM-scan groups along the 3rd dimension).
%
% Input
% -----
% cplxData:
%     interleaved, laterally-cropped complex volume for one channel
%     (depth x A-scans x B-scans), grouped in blocks of numBMscans
%     repeated B-scans along the 3rd dimension.
% numBMscans:
%     number of repeated B-scans per group (e.g. 4).
% usfac:
%     up-sampling factor for mcorrLocal sub-pixel registration.
%
% Output
% -----
% avgOCT:
%     mean B-scan amplitude per group.
% subOCTA:
%     sum of magnitudes of successive bulk-phase-corrected complex
%     differences per group.

numPoints = size(cplxData,1);
numAscans = size(cplxData,2);
numBscans = size(cplxData,3);

%% Sub-pixel Motion Correction %%

OCT_lmcorr = zeros(numPoints, numAscans, numBscans, 'like', cplxData);

for I = 1:numBMscans:numBscans
    OCT_lmcorr(:,:,I:I+numBMscans-1) = mcorrLocal(cplxData(:,:,I:I+numBMscans-1), usfac);
end

clear cplxData

%% OCTA Calculation %%

numGroups = numBscans / numBMscans;
avgOCT  = zeros(numPoints, numAscans, numGroups);
subOCTA = zeros(numPoints, numAscans, numGroups);

for I = 1:numBMscans:numBscans

    K = ((I-1)/numBMscans) + 1;

    BMscan = OCT_lmcorr(:,:,I:I+numBMscans-1);

    % Bulk phase correction relative to first B-scan
    Xconj   = BMscan(:,:,2:end) .* conj(BMscan(:,:,1));
    BulkOff = angle(sum(Xconj,1));
    BMscan(:,:,2:end) = BMscan(:,:,2:end) .* exp(-1j*BulkOff);

    % Average OCT
    avgOCT(:,:,K) = mean(abs(BMscan), 3);

    % Subtraction
    Diff = diff(BMscan, 1, 3);
    subOCTA(:,:,K) = sum(abs(Diff), 3);

end

end
