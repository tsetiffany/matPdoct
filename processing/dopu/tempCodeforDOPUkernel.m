%% Get the Thickness data from DOPU generation %%

Mean_table                      = []; %% Mean value of A lines
Mean_Rigid                      = []; %% Register all A line of single volume

for k = 1 : 33  %number of volumes
    Rigid_files                 = dir('*Thickness.mat');
    Rigid_dats                  = {Rigid_files.name}';
    Rigid_values                = load(Rigid_dats{k});
    for J = 1 : 500 %num alines per bscan
        Mean_Aline              = mean(Rigid_values.lowThickness(:, J));
        Mean_Rigid              = [Mean_Rigid , Mean_Aline];
    end
%     Mean_Rigid                  = Mean_Rigid';
    Mean_table                  = [Mean_table; Mean_Rigid , convertCharsToStrings(Rigid_dats{k})];
%     Mean_table                  = [Mean_Rigid , convertCharsToStrings(Rigid_dats{k}) ; Mean_table];
    Mean_Rigid                  = [];
    
end

% A is the DOPU values of Alines in single volume
% B is the name tag indicating which file contains the value
A                               = [];
B                               = [];
A                               = Mean_table(: , 1 : 500); %% Take the mean value of low DOPU regions
A                               = A';
B                               = Mean_table(: , 501); %% 
B                               = B';
Mean_Result                     = [B ; A];
% Mean Result giving you string format 

%%


A = Mean_table(: , 1 : 500); % Take the value and change it into double
A = A';

a = str2double(Mean_Result);
b = str2double(Mean_Result);

% Save the data to excel sheet (was tried to plot on excel)
% xlswrite('Adaptive_with10chunk_ker35.xls' , Mean_Result);
% xlswrite('Rigid.xls' , Mean_Result);



%%
% AT              = readtable('Adaptive_Cummulated_Alines.xls'); % Load the excel sheet data I don't think you need this
% RT              = readtable('Rigid_Cummulated_Alines.xls');
% T               = readtable('Adaptive_with10chunk_Cummulated_Alines.xls');
% T               = readtable('Adaptive_with10chunk_ker35_Cummulated_Alines.xls');
AT              = A;
ArrayTable      = table2array(AT);
ArrayTable      = ArrayTable(2 : 2000 , :); %four volumes in same directionality
ArrayTable_RG   = table2array(RT);
ArrayTable_RG   = ArrayTable_RG(2 : 2000 , :);



% Measure the mean and standard deviation values of all Alines

% Std_Array       = zeros(1,9);
% Mean_Array      = zeros(1,9);
% Std_Array(2 : 8)       = std(ArrayTable);
% Mean_Array(2 : 8)      = mean(ArrayTable , 1);
Std_Array       = std(ArrayTable);
Mean_Array       = mean(ArrayTable);
% Std_Array_Chunk = std(ArrayTable);
% Mean_Array_Chunk= mean(ArrayTable , 1);


% Std_Array_RG    = zeros(1,9);
% Mean_Array_RG   = zeros(1,9);
% Std_Array_RG(2 : 8)    = std(ArrayTable_RG);
% Mean_Array_RG(2 : 8)   = mean(ArrayTable_RG , 1);
Std_Array_RG       = std(ArrayTable_RG);
Mean_Array_RG       = mean(ArrayTable_RG , 1);

%% mean error bar plot


% Make X axis for directionality
% Plot errorbar with mean and standard deviation
x = [35, 25 , 15 , 5 , 0 , -5 , -15 , -25 , -35];
x = flip(x);
errorbar(x , Mean_Array_RG , Std_Array_RG);
hold on; 

errorbar(x , Mean_Array , Std_Array , 'r');

err = 8*ones(size(y));
errorbar(x,y,err)


%% Max and min plot 
MaxRV = [];
MinRV = [];
MaxAV = [];
MinAV = [];


for K = 1 : 7
    
    MaxAV = [MaxAV , max(ArrayTable(: , K))];
    MinAV = [MinAV , min(ArrayTable(: , K))];
    MaxRV = [MaxRV , max(ArrayTable_RG(: , K))];
    MinRV = [MinRV , min(ArrayTable_RG(: , K))];
    
end
x = [25 , 15 , 5 , 0 , -5 , -15 , -25];
x = flip(x);

errorbar(x,(MinRV + MaxRV)/2,(MaxRV - MinRV)/2);

hold on
errorbar(x , (MinAV + MaxAV)/2 , (MaxAV - MinAV)/2 , 'r');

%%


y = rand(1,10); % your mean vector;

x = 1:numel(y);

std_dev = 1;

curve1 = y + Std_Array_RG;

curve2 = y - Std_Array_RG;

x2 = [x, fliplr(x)];

inBetween = [curve1, fliplr(curve2)];

fill(x2, inBetween, 'b');

hold on;

plot(x, y, 'r', 'LineWidth', 2);