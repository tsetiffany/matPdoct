function u12u12 = unpack_u12u16(filename,numPoints,numAscans, ref_Bscan)

total_elem = numPoints * numAscans;
byteSize = 12 / 8;

fid = fopen(filename,'r');
fseek(fid,byteSize*numPoints*numAscans*(ref_Bscan-1),-1);
data = fread(fid,total_elem*byteSize,"uint8=>int32");
fclose(fid);

u12u12 = zeros(total_elem,1,'int32');
u12u12(1:2:end) = data(1:3:end);
u12u12(2:2:end) = bitshift(data(3:3:end),4);

u12u12(1:2:end) = u12u12(1:2:end) + bitshift(bitand(data(2:3:end), hex2dec('F')), 8);
u12u12(2:2:end) = u12u12(2:2:end) + bitshift(bitand(data(2:3:end), hex2dec('F0')), -4);

u12u12 = bitshift(u12u12,4);

u12u12 = reshape(u12u12,numPoints,numAscans);

% plot(u12u12(:,1,1))