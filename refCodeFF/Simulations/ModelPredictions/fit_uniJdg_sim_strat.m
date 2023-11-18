function estimatedP = fit_uniJdg_sim_strat(binary_unityJdg, s_A, s_V, a_A,...
    b_A, sigma_AV_A, sigma_AV_V, CI, mu_P, lapse, model)
global lenC lenP
switch model.strategy_unity
    case 'PosteriorC1'
        lb = ones(1,lenC*lenP).*0.01; plb = ones(1,lenC*lenP).*0.1;
        ub = ones(1,lenC*lenP).*0.99; pub = ones(1,lenC*lenP).*0.9;
        init = ones(1,lenC*lenP).*0.5;
        nLogL = @(p) sum(arrayfun(@(idx) get_mPred(p(idx), NaN, s_A, s_V, a_A,...
            b_A,sigma_AV_A,sigma_AV_V, CI, mu_P, lapse, model,...
            squeeze(binary_unityJdg(ceil(idx/2),rem(idx+1,2)+1,:,:,:))),1:lenC*lenP));
        [estimatedP,~] = bads(nLogL, init, lb, ub, plb,pub);
    otherwise
        lb = ones(1,lenC*lenP+1).*0.01;      plb = [ones(1,lenC*lenP).*0.1,3];
        ub = [ones(1,lenC*lenP).*0.99,15]; pub = [ones(1,lenC*lenP).*0.9,10];
        init = [ones(1,lenC*lenP).*0.9,5];
        nLogL = @(p) sum(arrayfun(@(idx) get_mPred(p(idx), p(lenC*lenP+1), s_A, s_V, a_A,...
            b_A,sigma_AV_A,sigma_AV_V, CI, mu_P, lapse, model,...
            squeeze(binary_unityJdg(ceil(idx/2),rem(idx+1,2)+1,:,:,:))),1:lenC*lenP));
        [estimatedP,~] = bads(nLogL, init, lb, ub, plb,pub);
end
        
        
        