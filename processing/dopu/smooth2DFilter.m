function result = smooth2DFilter(input, kernel)

imgDepth = size(input,1);
imgWidth = size(input,2);

weight_Width_S = linspace(1,kernel(2)-1,kernel(2)-1)';
weight_Width_M = ones([imgWidth-2*(kernel(2)-1), 1])*kernel(2);
weight_Width_E = flipud(weight_Width_S);
weight_Width   = cat(1,weight_Width_S,weight_Width_M,weight_Width_E);
% e.g., kernel [1,8], input 372 pts by 1000 alines, weight_width =
% [1,2,3...8....8....3,2,1] for 1000 values.  Convolution kernels are all
% ones adn then weighted by the weight_width afterwards
    
for k = 1:imgDepth
    kernelY  = ones([1 kernel(2)]);
    img_conY = conv(input(k,:), kernelY, 'same');
    outputY(k,:) = img_conY'./weight_Width;
end

weight_Depth_S = linspace(1,kernel(1)-1,kernel(1)-1)';
weight_Depth_M = ones([imgDepth-2*(kernel(1)-1), 1])*kernel(1);
weight_Depth_E = flipud(weight_Depth_S);
weight_Depth   = cat(1,weight_Depth_S,weight_Depth_M,weight_Depth_E);

outputX = outputY;
for kk = 1:imgWidth
    kernelX  = ones([kernel(1),1]);
    img_conX = conv(outputY(:,kk), kernelX, 'same');
    outputX(:,kk) = img_conX./weight_Depth;
end

result = outputX;

end