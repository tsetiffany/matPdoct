function imout = imcenhance(img)

% imout = adapthisteq(mat2gray(img));
imout = (mat2gray(img));
s = size(img,1);
imout = imresize(imout, [2*s 2*s]);
imout = notchfilter(imout);
imout = imadjust(mat2gray(imout));
% imshow(imout)

end