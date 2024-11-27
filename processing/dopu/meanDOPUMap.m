function [RPE, choroid] = meanDOPUMap(DOPU,DOPU_depth_C, DOPU_depth_C_bot, depth)
   % Takes the filtered, upright DOPU volume and the DOPU depth lines (top
   % boundary and bottom boundary) and produces mean maps for the RPE (with
   % thickness of "depth" pixels, usually about 5 pixels) and choroid (the rest)
   
    % size values
    numPoints=size(DOPU_test,1);
    numAlines=size(DOPU_test,2);
    numBscans=size(DOPU_test,3);
   
    % prepare projections
    RPE = zeros(numAlines, numBscans);
    choroid = zeros(numAlines, numBscans);
   
    for B=1:numBscans % loop across the frames
        for A=1:numAlines % loop across the A-lines
            % extract the indices (and round them)
            top = round(DOPU_depth_C(A,B));
            mid = top+depth;
            bot = round(DOPU_depth_C_bot(A,B));
            
            % ensure no out of bounds conditions
            if mid >numPoints
                mid=numPoints;
            end
            
            % extract the mean values
            RPE(A,B) = mean(DOPU(top:mid,A,B));
            choroid(A,B) = mean(DOPU(mid:bot,A,B);
        end
    end


end