function testRect(type)
% ----------------------------------------------------------------------
% testRect(type)
% ----------------------------------------------------------------------
% Goal of the function :
% Simple function that illustrates the FillRect and FrameRect sub
% functions of Screen()
% ----------------------------------------------------------------------
% Input(s) :
% type : switch of the display
% if type == 1 : FILLRECT
% if type == 2 : FRAMERECT
% else : FILLRECT & FRAMERECT
% ----------------------------------------------------------------------
% Output(s):
% (none)
% ----------------------------------------------------------------------
% Function created by Martin SZINTE (martin.szinte@gmail.com)
% Last update : 02 / 01 / 2014
% Project : Programming courses
% Version : -
% ----------------------------------------------------------------------
type;
scrAll = Screen('Screens');
scrNum = max(scrAll);
HideCursor;
Screen('Preference', 'SkipSyncTests', 1);
[scr] = Screen('OpenWindow',scrNum);
priorityLevel = MaxPriority(scr);Priority(priorityLevel);
colBG = [255 255 255]; % screen background color
colRect = [0 0 0]; % filling color of the rectangles
rectRect = [200, 200, 400, 300]; % [Left Top Right Bottom] coordinates of the first
rectangle
rectRect2 = [200, 400, 700, 500]; % [Left Top Right Bottom] coordinates of the second
rectangle
rectWidth = 20; % width of the frame rectangle

% Main display loop
for timeFlip = 1:200

 Screen('FillRect',scr,colBG); % fills background color of the screen

 % display selection loop
 if type == 1
 Screen('FillRect',scr,colRect,rectRect); % draws a filled rectangle on the secondary screen buffer
 elseif type == 2
 Screen('FrameRect',scr,colRect,rectRect2,rectWidth); % draws a framed rectangle on the secondary screen buffer
 else
 Screen('FillRect',scr,colRect,rectRect); % draws a filled rectangle on the secondary screen buffer
 Screen('FrameRect',scr,colRect,rectRect2,rectWidth); % draws a framed rectangle on the secondary screen buffer
 end

 Screen(scr, 'Flip'); % Flips the secondary screen buffer on the monitor

end
ShowCursor;
Screen('CloseAll');
end