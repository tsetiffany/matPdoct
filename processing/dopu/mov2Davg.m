function output = mov2Davg(input, window)

for i = 1:size(input,3)
    output(:,:,i) = imguidedfilter(input(:,:,i),'NeighborhoodSize', window);
end
