function estimatedP = fit_uniJdg_sim(binary_unityJdg, s_A, s_V, a_A,...
    b_A, sigma_AV_A, sigma_AV_V, epsilon,CI, mu_P, lapse, model, MType)
global lenC lenP
estimatedP = NaN(lenC, lenP);
p = 1e-2:1e-1:0.99;

switch MType
    case 'full'
        for i = 1:lenC
            for j = 1:lenP
                D = squeeze(binary_unityJdg(i,j,:,:,:));
                nLogL = @(p) get_mPred(p,epsilon, s_A, s_V, a_A, b_A,sigma_AV_A,...
                    sigma_AV_V, CI, mu_P, lapse, model, D);
                nLogL_allp = arrayfun(@(idx) nLogL(p(idx)), 1:length(p));
                [~, min_idx] = min(nLogL_allp);
                estimatedP(i,j) = p(min_idx);
            end
        end
    case 'restricted'
        D_pre1 = squeeze(binary_unityJdg(1,1,:,:,:));
        D_pre2 = squeeze(binary_unityJdg(2,1,:,:,:));
        nLogL_pre = @(p) get_mPred(p,epsilon, s_A, s_V, a_A, b_A,sigma_AV_A, ...
            sigma_AV_V, CI, mu_P, lapse, model, D_pre1) + ...
            get_mPred(p, epsilon,s_A, s_V, a_A, b_A, sigma_AV_A, sigma_AV_V,...
            CI, mu_P, lapse, model, D_pre2);
        nLogL_allpre = arrayfun(@(idx) nLogL_pre(p(idx)), 1:length(p));
        [~, min_idx] = min(nLogL_allpre);
        estimatedP(:,1) = p(min_idx);
        for i = 1:lenC
            D_post = squeeze(binary_unityJdg(i,2,:,:,:));
            nLogL_post = @(p) get_mPred(p, epsilon,s_A, s_V, a_A, b_A, sigma_AV_A,...
                sigma_AV_V, CI, mu_P, lapse, model, D_post);
            nLogL_post = arrayfun(@(idx) nLogL_post(p(idx)), 1:length(p));
            [~, min_idx] = min(nLogL_post);
            estimatedP(i,2) = p(min_idx);
        end      
    case 'null'
        for i = 1:lenC
            D_phase1 = squeeze(binary_unityJdg(i,1,:,:,:));
            D_phase2 = squeeze(binary_unityJdg(i,2,:,:,:));
            nLogL_cond = @(p) get_mPred(p,epsilon, s_A, s_V, a_A, b_A, sigma_AV_A,...
                sigma_AV_V, CI, mu_P, lapse, model, D_phase1) + ...
                get_mPred(p,epsilon, s_A, s_V, a_A, b_A,sigma_AV_A, sigma_AV_V,...
                CI, mu_P, lapse, model, D_phase2);
            nLogL_allphases = arrayfun(@(idx) nLogL_cond(p(idx)), 1:length(p));
            [~, min_idx] = min(nLogL_allphases);
            estimatedP(i,:) = p(min_idx);            
        end
end
end
        
        
        