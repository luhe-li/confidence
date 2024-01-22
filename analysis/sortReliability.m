clear;

load('uniLoc_sub2_ses-V.mat');

[~, temp1] = sort(VSinfo.SD_blob);
reliSortResp = Resp(temp1);
reli1resp = reliSortResp(1:ExpInfo.nRep * ExpInfo.nLevel);
reli2resp = reliSortResp((ExpInfo.nRep * ExpInfo.nLevel+1):end);

[~, temp2] = sort([reli1resp.target_idx]);
sortedReli1Resp = reli1resp(temp2);
[~, temp3] = sort([reli2resp.target_idx]);
sortedReli2Resp = reli2resp(temp3);

save('uniLoc_sub2_ses-V2.mat','Resp','reliSortResp','sortedResp','ExpInfo','ScreenInfo','VSinfo','AudInfo','sortedReli1Resp','sortedReli2Resp')