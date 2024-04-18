clear;
sA = 1:10;
sV = -10:1:-1;
sAV = combvec(sA, sV);
numS = numel(sAV);


x = 1:100000;
xr = reshape(x, repmat(10, 1, 5));
