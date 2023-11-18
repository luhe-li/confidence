function sim_pC1 = sim_reportingC1(pReportingC1, nT_perPair)

global lenC lenP lenS
randN = rand([lenC, lenP, lenS, lenS, nT_perPair]);
pReportingC1 = repmat(pReportingC1, [1,1,1,1,nT_perPair]);
bool_less = randN < pReportingC1;
sim_pC1 = sum(bool_less,5)./nT_perPair;


