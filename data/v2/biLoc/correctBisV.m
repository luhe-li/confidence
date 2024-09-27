clear;

for sub = 11:14

flnm = sprintf('biLoc_sub%i_ses1', sub);

load(flnm);

% update visual locations in different units
ExpInfo.randVisCM     = ExpInfo.randVisPixel ./ ScreenInfo.numPixels_perCM;
ExpInfo.randVisVA    = (180/pi) * (atan(ExpInfo.randVisCM ./ExpInfo.sittingDistance));

save(flnm)

end