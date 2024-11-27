function parameters = getParameters(fn)
fn = fn(1:end-4);
fn_xml = strcat(fn,'.xml');
file_xml = fileread(fn_xml);
file_xml(strfind(file_xml, '=')) = [] ;  % Omit all ':' character from text file
file_xml(strfind(file_xml, '"')) = [] ;  % Omit all ':' character from text file
key1    = 'Width';
key2    = 'Height';
key3    = 'Number_of_Frames';
key4    = 'Number_of_Volumes';
key5    = 'Number_of_BM_scans';
key6    = 'C2';
key7    = 'C3';

Idx1    = strfind(file_xml, key1);
Idx2    = strfind(file_xml, key2);
Idx3    = strfind(file_xml, key3);
Idx4    = strfind(file_xml, key4);
Idx5    = strfind(file_xml, key5);
Idx6    = strfind(file_xml, key6);
Idx7    = strfind(file_xml, key7);


numPoints    = sscanf(file_xml(Idx1(1)+length(key1):end), '%g', 1);
numAscans    = sscanf(file_xml(Idx2(1)+length(key2):end), '%g', 1);
numBscans    = sscanf(file_xml(Idx3(1)+length(key3):end), '%g', 1);
numCscans    = sscanf(file_xml(Idx4(1)+length(key4):end), '%g', 1);
numMscans    = sscanf(file_xml(Idx5(1)+length(key5):end), '%g', 1);
disCoeff2    = sscanf(file_xml(Idx6(1)+length(key6):end), '%g', 1);
disCoeff3    = sscanf(file_xml(Idx7(1)+length(key7):end), '%g', 1);


parameters = [numPoints numAscans numBscans numCscans numMscans disCoeff2 disCoeff3];

end
