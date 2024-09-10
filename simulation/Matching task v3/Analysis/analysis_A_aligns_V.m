%% This script analyzes the AV alignment data
clear all;close all;clc
%-------------------------------------------------------------------------
%load the data
%-------------------------------------------------------------------------
subjNum              = 21;
subjInitial          = 'LL';
C                    = load(strcat('A_aligns_V_sub', num2str(subjNum), '.mat'),...
                        'A_aligns_V_data');
ExpInfo              = C.A_aligns_V_data{1};
%get the updating distance for all the staircases and convert it to deg
Distance             = C.A_aligns_V_data{8};
Distance(:,end)      = []; %the last one was generated but never used in the experiment
%get the response (whether participants think the V is to the left or right of the A)
LeftOrRight          = C.A_aligns_V_data{6};
%get the data for easy trials
data_easyTrials      = C.A_aligns_V_data{10};

%-------------------------------------------------------------------------
%Get relevant information
%-------------------------------------------------------------------------

%there are 4 locations of the V: [-7.5, -2.5, 2.5, 7.5]
locations_V         = C.A_aligns_V_data{3}.locations_deg;
%for each test location, there were 2 staircases
%1. the A starts from the left side of the V
%2. the A starts from the right side of the V
numStaircasesPerLoc = ExpInfo.numStaircases/ExpInfo.testLocations;
%number of easy trials
numEasyTrials       = length(data_easyTrials);
%number of easy trials for each test location
numEasyTrialsPerLoc = size(data_easyTrials,2)/length(locations_V);

%% plot the interleaved staircases
figure(1)
%create a color matrix for those 8 interleaved staircases
colorMat = [70,130,180;255,140,0; 34,139,34;255,0,0]./255;
%for each test location, there are 2 interleaved staircases
for k = 1:ExpInfo.testLocations
    for i = 1:numStaircasesPerLoc
        %plot the location of the A
        plot(1:ExpInfo.numTrials,Distance((k-1)*2+i,:),'-','LineWidth',2,...
            'Color',colorMat(k,:)); hold on
    end
    %plot the location of the V
    plot([1 ExpInfo.numTrials],[locations_V(k) locations_V(k)],'--',...
        'LineWidth',2,'Color',colorMat(k,:)); hold on
    
    %plot the responses (blue: A is to the left of V; yellow: A is to the right of V)
    for i = 1:numStaircasesPerLoc
        for j = 1:ExpInfo.numTrials 
            if LeftOrRight((k-1)*2+i,j)==-1
                plot(j,Distance((k-1)*2+i,j),'Marker','s','MarkerEdgeColor',...
                    'y','MarkerFaceColor','b','MarkerSize',7);
            else
                plot(j,Distance((k-1)*2+i,j),'Marker','s','MarkerEdgeColor',...
                    'b','MarkerFaceColor','y','MarkerSize',7);
            end
            hold on
        end
    end
    xlabel('Trial number'); ylabel('Stimulus location (deg)');
end
hold off; box off
xlim([1 ExpInfo.numTrials]);
%xticks([1 ExpInfo.numTrials/2 ExpInfo.numTrials]); yticks(testLocs);
title('Interleaved Staircases');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.7]);
set(gca,'FontSize',15, 'XTick',[1 ExpInfo.numTrials/2 ExpInfo.numTrials],...
    'YTick',locations_V);
saveas(gcf,sprintf(['InterleavedStaircases_sub' num2str(subjNum) '.jpg']));

%% sort the data and add easy trials
%combine the loc of the V and the responses for the 2 interleaved staircases 
%for each test location together
Distance_reshaped    = reshape(Distance',[ExpInfo.numTrials*numStaircasesPerLoc,...
                        ExpInfo.testLocations])';
LeftOrRight_reshaped = reshape(LeftOrRight',[ExpInfo.numTrials*numStaircasesPerLoc,...
                        ExpInfo.testLocations])';
numTrials_comb       = size(Distance_reshaped,2);

%organize the data
for i = 1:ExpInfo.testLocations
    %put the location of the V in the 1st row and the response in the 2nd row
    data.sorted{i}(1,1:numTrials_comb) = Distance_reshaped(i,:);
    data.sorted{i}(2,1:numTrials_comb) = LeftOrRight_reshaped(i,:);
    
    %from easy trials, find the indices of the trials that match the test location
    matching_idx = (data_easyTrials(3,:) == locations_V(i));
    %append the location of the A to data.sorted
    data.sorted{i}(1,numTrials_comb+1:numTrials_comb+numEasyTrialsPerLoc) = ...
        data_easyTrials(4,matching_idx);
    %append the response to data.sorted
    data.sorted{i}(2,numTrials_comb+1:numTrials_comb+numEasyTrialsPerLoc) = ...
        data_easyTrials(6,matching_idx);
   
    %when fitting a psychometric function, we will calculate the
    %probability of saying the V is to the right of the A, so we change all
    %-1 to 0
    data.sorted{i}(2,data.sorted{i}(2,:)==-1) = 0;
end

%% fit a psychometric function
%common lapse rate, alpha1, sigma1, alpha2, sigma2, alpha3, sigma3, alpha4,sigma4
lb              = [   0, -10, 0.01,  -4, 0.01,  2, 0.01,  7, 0.01]; 
ub              = [0.06,   0,    5,   4,    6,  10,    3,  14,    3]; %set boundaries
initialization  = [   0,  -5,    1,  -1,    1,   7,    1,  10,    1];
options         = optimoptions(@fmincon,'MaxIterations',1e5,'Display','off');
nLogL           = @(p) nLL_commonLapse(p(1),p(2),p(3),p(4),p(5),p(6),p(7),...
                    p(8),p(9),data.sorted);
[data.estimatedP, LogLValue] = fmincon(nLogL, initialization,[],[],[],[],...
                    lb,ub,[],options);
disp('Estimated values of the parameters:');
disp(data.estimatedP);

%% plot fitted psychometric function 
%specify plotting information
bds               = [min(Distance(:)), max(Distance(:))]; numX = 1e3;
pltInfo.x         = linspace(bds(1)-3,bds(end)+3, numX);
pltInfo.numBins   = 30; 
pltInfo.bool_save = 0;
pltInfo.subjI     = subjInitial;
data.binned       = plot_psychfunc(data, ExpInfo, pltInfo);

%% Do bootstrap and calculate the confidence interval of PSE
numBootstraps                      = 1e3;
[data.PSE,data.PSE_lb,data.PSE_ub] = Bootstrap(data.sorted,numBootstraps,...
                                        lb,ub,options);
disp(mean(data.PSE));

%% plot the PSEs 
%fit the data with a line
data.polyfit      = polyfit(locations_V,data.estimatedP(2:2:end),1);
%calculate the predicted value
data.polyval      = polyval(data.polyfit,locations_V);
%calculate the residuals
y_residuals       = data.estimatedP(2:2:end) - data.polyval;
%calculate the sum squared error
Sum_squaredErrors = sum(y_residuals.^2);
Sum_squaresTotal  = sum((data.estimatedP(2:2:end) - mean(data.estimatedP(2:2:end))).^2);
%calculate R square
data.R_square     = 1 - Sum_squaredErrors/Sum_squaresTotal;
disp(['R^2: ' num2str(round(data.R_square,4))]);
%calculate RMSE
data.RMSE         = sqrt(Sum_squaredErrors/(length(data.sorted{i})*length(locations_V)-1));
disp(['RMSE: ' num2str(round(data.RMSE,4))]);
%plot it
plot_PSE(data, ExpInfo, pltInfo)

%% save the data
out1FileName = ['AV_alignment_sub' num2str(subjNum) '_dataSummary'];
AV_alignment_data = {ExpInfo, data};
save(out1FileName,'AV_alignment_data');
