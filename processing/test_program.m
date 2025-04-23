
function test_program(msize)

%generate a matrix of random numbers of dimension msize x msize

rmatrix=rand(msize);

%create an output file name based on msize and write the matrix to it

fname=num2str(msize);

fname=strcat(fname,'.csv');

writematrix(rmatrix, fname);

end