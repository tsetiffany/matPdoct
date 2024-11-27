function PostprocScript(folder)

folder = 'C:\Users\MJJU\Desktop\BiDiOCTA';


addpath('D:\Dropbox\ProgramScripts\MatlabScripts\AO-OCTA\AO-OCTA_ver3');
mcorr_path = strrep(folder,'RAW','mcorr');

rescaleFolder = 'C:\Users\MJJU\Desktop\BiDiOCTA';
fn_ResParam = 'LUTSS.dat';
fid_ResParam = fopen(fullfile(rescaleFolder,fn_ResParam));
rescaleParam = fread(fid_ResParam, 'double')+1.00;
rescaleParam =  rescaleParam;
fclose(fid_ResParam);


options.rescaleParam    = rescaleParam;
options.dispMaxOrder    = 5;
options.coeffRange      = 100;

%%%%% Find filenames of RAW files to be processed %%%%%
cd(folder);
files   = (dir('*.unp'));
fnames  = {files.name}';

%%%%% Cplx OCT processing start %%%%%
for K = 1:size(fnames,1)
    proCplxOCT(fnames{K}, options,folder);
    proCplxOCA(fnames{K}, options);
end
    