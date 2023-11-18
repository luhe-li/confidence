%In this task, participants practiced using the rollerball mouse to adjust 
%the horizontal location of a response cursor to localize a stimulus. Each 
%trial started with a central fixation cross presented for 500 ms. Then, 
%the visual stimulus, a small white square, was displayed on the screen for 
%100 ms, followed by a blank screen presented for 1000 ms. Next, a response 
%cursor appeared on the screen. Participants moved the response cursor to 
%the location of the stimulus by moving the rollerball, and clicked the 
%mouse to register the response. There was no time limit for the response. 
%Visual feedback of the cursor location was provided during adjustment, but 
%error feedback was not provided (Fig. 3A). There were eight possible 
%horizontal positions for the stimulus, evenly spaced from -17.5 to 17.5◦ 
%in steps of 5◦. Each stimulus location was visited 30 times in random 
%order, resulting in a total of 240 trials. The inter-trial interval was 
%500 ms. This experiment usually took half an hour to complete.

%This is a practice

%% Enter subject's name
clear all; close all; clc
rng('shuffle');

ExpInfo.subjID = [];
while isempty(ExpInfo.subjID) == 1
    try ExpInfo.subjID = input('Please enter participant ID#: '); %'s'
    catch
    end
end
out1FileName = ['PointingTest_practice_sub' num2str(ExpInfo.subjID)];

%% Screen Setup 
Screen('Preference', 'VisualDebugLevel', 1);
Screen('Preference', 'SkipSyncTests', 1);
[windowPtr,rect] = Screen('OpenWindow', 0, [0,0,0]);
%[windowPtr,rect] = Screen('OpenWindow', 0, [105,105,105],[100 100 1000 480]); % for testing
[ScreenInfo.xaxis, ScreenInfo.yaxis] = Screen('WindowSize',windowPtr);
Screen('TextSize', windowPtr, 35) ;   
Screen('TextFont',windowPtr,'Times');
Screen('TextStyle',windowPtr,1); 

[center(1), center(2)]     = RectCenter(rect);
ScreenInfo.xmid            = center(1); % horizontal center
ScreenInfo.ymid            = center(2); % vertical center
ScreenInfo.backgroundColor = 105;
ScreenInfo.numPixels_perCM = 7.5;
ScreenInfo.liftingYaxis    = 300; 

%% specify the experiment informations
ExpInfo.stimulusLocs_deg = -17.5:5:17.5;
ExpInfo.sittingDistance  = 105; %the distance between the screen and participants
ExpInfo.stimulusLocs_cm  = round(tan(deg2rad(ExpInfo.stimulusLocs_deg)).*...
                            ExpInfo.sittingDistance,1); %in cm 
ExpInfo.numTrialsPerLoc  = 1;
ExpInfo.numTrials        = ExpInfo.numTrialsPerLoc * length(ExpInfo.stimulusLocs_deg); 
ExpInfo.numFrames_target = 6; %100ms

%date_easyTrials stores all the data for randomly inserted easy trials
%1st row: the target will appear at one of those 10 locations (pseudorandomized)
%2nd row: response (degree)
%3rd row: the tarlget location (in pixel)
%4th row: response (in pixel)
%5th row: response time
data      = zeros(2,ExpInfo.numTrials);
data(1,:) = Shuffle(Shuffle(repmat(ExpInfo.stimulusLocs_deg,...
                [1 ExpInfo.numTrialsPerLoc])));
data(3,:) = round(ScreenInfo.numPixels_perCM.*tan(deg2rad(data(1,:))).*...
                ExpInfo.sittingDistance + ScreenInfo.xmid);

%% Start the experiment
%record movie
% moviePtr = Screen('CreateMovie',windowPtr,'PointingPractice.mov',...
%     ScreenInfo.xaxis,ScreenInfo.yaxis,60);
DrawFormattedText(windowPtr, 'Press any button to start the experiment.',...
    'center',ScreenInfo.yaxis-ScreenInfo.liftingYaxis,[255 255 255]);
Screen('Flip',windowPtr);
%Screen('AddFrameToMovie',windowPtr);
KbWait(-3); %-1 if using the keyboard; -3 if using the attached keypad
WaitSecs(1);
Screen('Flip',windowPtr);

for i = 1:ExpInfo.numTrials 
    Screen('FillRect', windowPtr,[255 255 255], [ScreenInfo.xmid-7 ...
        ScreenInfo.yaxis-ScreenInfo.liftingYaxis-1 ScreenInfo.xmid+7 ...
        ScreenInfo.yaxis-ScreenInfo.liftingYaxis+1]);
    Screen('FillRect', windowPtr,[255 255 255], [ScreenInfo.xmid-1 ...
        ScreenInfo.yaxis-ScreenInfo.liftingYaxis-7 ScreenInfo.xmid+1 ...
        ScreenInfo.yaxis-ScreenInfo.liftingYaxis+7]);
    %display how much have subjects done in percentage
    DrawFormattedText(windowPtr,[num2str(i) '/' num2str(ExpInfo.numTrials)],...
        'center', ScreenInfo.yaxis-ScreenInfo.liftingYaxis-50,...
        [255 255 255]);
    Screen('Flip',windowPtr);
    WaitSecs(1);
    Screen('Flip',windowPtr);
    WaitSecs(0.5);
    
    [data(2,i), data(4,i), data(5,:)] = PresentVisualStimulus(data(3,i),...
        ScreenInfo, ExpInfo, windowPtr);
end
ShowCursor;
%Screen('FinalizeMovie',moviePtr);

%% Save data and end the experiment
PointingTest_practice_data = {ExpInfo,ScreenInfo, data};
save(out1FileName,'PointingTest_practice_data');
Screen('CloseAll');
