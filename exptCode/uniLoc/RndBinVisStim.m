function Resp = RndBinVisStim(i, ExpInfo, ScreenInfo,VSinfo,windowPtr)

%% precompute visual stimuli

%first compute the location of the visual stimulus in pixels
loc_pixel = round(ExpInfo.randVisPixel(i));
targetLoc = [ScreenInfo.xmid + loc_pixel,300];
num_dots = 1000;
RNcoordinates = round(rand(num_dots,2) .* repmat([1024,768],num_dots,1));

coherPDF = mvnpdf(RNcoordinates,targetLoc,[400,0;0,400]);
coherPDF = coherPDF ./ sum(coherPDF) .*160;
num_frames = 180;
frameBin = randi([0 1],num_dots,num_frames);
ideal = repmat([0;0;1;1],num_frames./4,1);
include = 0;
for f = 1:num_frames
    randInx = rand(num_dots,1)<coherPDF;
    frameBin(randInx,f) = ideal(f);
    frameBin(~randInx,f) = randi([0 1],sum(~randInx),1);
    include = include + sum(randInx);
end
sort(coherPDF)
%%
filename = '0011.gif';

set(gcf, 'Color', [0.5 0.5 0.5]);  % Figure background color
set(gca, 'Color', [0.5 0.5 0.5]);  % Figure background color

h = figure(1);
colormap("gray");

for f = 1:num_frames
    clf;
    hold on
    rectangle('Position', [0, 0, 1024, 768], ...
          'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
    % plot(targetLoc(1),targetLoc(2),'ro')
    scatter(RNcoordinates(:,1),RNcoordinates(:,2),[],frameBin(:,f),'filled')
    xlim([0,1024])
    ylim([0,768])
    hold off
    % Capture the plot as an image
    frame = getframe(h);
    img = frame2im(frame);
    [imind, cm] = rgb2ind(img, 256);
    
    % Write to the GIF File
    if f == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.067);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.067);
    end
end

%% start the trial

% fixation
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('FillRect', windowPtr,[255 255 255], [ScreenInfo.x1_lb,...
    ScreenInfo.y1_lb, ScreenInfo.x1_ub, ScreenInfo.y1_ub]);
Screen('FillRect', windowPtr,[255 255 255], [ScreenInfo.x2_lb,...
    ScreenInfo.y2_lb, ScreenInfo.x2_ub, ScreenInfo.y2_ub]);
Screen('Flip',windowPtr);
WaitSecs(ExpInfo.tFixation);

% blank screen 1
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);
WaitSecs(ExpInfo.tBlank1);

% display visual stimulus
Screen('DrawTexture',windowPtr,dotClouds_targetLoc,[],...
        [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
vbl = Screen('Flip',windowPtr);
Screen('Flip',windowPtr, vbl + (VSinfo.numFrames - 0.5) * ScreenInfo.ifi);

% mask
for jj = 1:VSinfo.numFramesMasker 
Screen('DrawTexture', windowPtr, VSinfo.gwn_texture(rem(jj,VSinfo.GWNnumFrames)+1),[],...
         [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);
end

Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);

%% response

% perceptual response
yLoc = ScreenInfo.yaxis-ScreenInfo.liftingYaxis;
SetMouse(randi(ScreenInfo.xaxis*2,1), yLoc*2, windowPtr);
HideCursor;
resp = 1;
tic;
stopRecorded = 0;
x = -1;
while resp
    cache = x;
    [x,~,~] = GetMouse(windowPtr);
    HideCursor;
    x = min(x, ScreenInfo.xmid*2);
    x = max(0,x);
    Screen('DrawTexture',windowPtr, VSinfo.grey_texture,[],...
        [0,0,ScreenInfo.xaxis, ScreenInfo.yaxis]);
    Screen('DrawLine', windowPtr, [255 255 255],x, yLoc-3, x, yLoc+3, 1);
    Screen('Flip', windowPtr);
    
    locdiff = abs(cache - x);
    if locdiff ~= 0
        stopRecorded = 0;
    elseif locdiff == 0 && ~stopRecorded
        mouseStopT = GetSecs();
        stopRecorded = 1;
    end
    
    % Check the keyboard
    [keyIsDown, startTime, keyCode] = KbCheck();
    if keyIsDown
        % Check if any of the specified keys are pressed
        if keyCode(KbName('a')) || keyCode(KbName('s')) || keyCode(KbName('d')) || keyCode(KbName('f'))
            [releaseTime, ~, ~] = KbReleaseWait();
            Resp.PressDuration = releaseTime - startTime;
            Resp.mouseStopDuration = startTime - mouseStopT;
            if keyCode(KbName('a'))
                conf = 1;
            elseif keyCode(KbName('s'))
                conf = 2;
            elseif keyCode(KbName('d'))
                conf = 3;
            elseif keyCode(KbName('f'))
                conf = 4;
            end
            resp = 0;
        end
        if keyCode(KbName('ESCAPE'))
            sca;
            ShowCursor;
            Screen('CloseAll');
            error('Escape!');
        end
    end
end
Resp.conf = conf;
Resp.RT1  = toc;
Resp.response_pixel = x;
Resp.response_cm    = (Resp.response_pixel -  ScreenInfo.xmid)/ScreenInfo.numPixels_perCM;
Resp.response_deg   = rad2deg(atan(Resp.response_cm/ExpInfo.sittingDistance));
HideCursor;

% ITI
Screen('DrawTexture',windowPtr,VSinfo.grey_texture,[],...
    [0,0,ScreenInfo.xaxis,ScreenInfo.yaxis]);
Screen('Flip',windowPtr);
WaitSecs(ExpInfo.ITI);

% Record target location
Resp.target_idx = ExpInfo.randVisIdx(i); % visual location that corresponds to speaker index
Resp.target_cm = ExpInfo.randVisCM(i);
Resp.target_pixel = ExpInfo.randVisPixel(i);
Resp.target_deg = ExpInfo.randVisVA(i);

end