function out = NLL_MA_gauss(freeParam, model, data)

switch model.mode

    case 'initiate'

        out.param_id                  = {'aA','bA','\sigma_{V}','\sigma_{A}','\sigma_{P}','\sigma_{C}','p_{common}'};
        out.num_param                 = length(out.param_id);

        % hard bounds, the range for LB, UB, larger than soft bounds
        paraH.aA     = [-3, 3];       % scale
        paraH.bA     = [-10, 10];      % intercept
        paraH.sigV   = [1e-2, 5];      % cm
        paraH.sigA   = [1, 20];        % cm
        paraH.sigP   = [3, 30];
        paraH.sigC   = [0.01, 10];      % measurement noise of confidence
        paraH.pC1    = [1e-3, 1-1e-3]; % weight

        % soft bounds, the range for PLB, PUB
        paraS.aA     = [-0.5, 1.5];      % scale
        paraS.bA     = [-5, 5];      % intercept
        paraS.sigV   = [1, 3];         % cm
        paraS.sigA   = [5, 10];        % cm
        paraS.sigP   = [5, 10];
        paraS.sigC   = [1, 5];         % measurement noise of confidence
        paraS.pC1    = [0.4, 0.6];     % weight

        % reorganize parameter bounds to feed to bads
        fn                           = fieldnames(paraH);
        for k                        = 1:numel(fn)
            out.lb(:,k)                  = paraH.(fn{k})(1);
            out.ub(:,k)                  = paraH.(fn{k})(2);
            out.plb(:,k)                 = paraS.(fn{k})(1);
            out.pub(:,k)                 = paraS.(fn{k})(2);
        end
        out.paraS                    = paraS; out.paraH = paraH;

        % get grid initializations
        model.num_sec = model.n_run*2;
        out.init                     = getInit(out.lb, out.ub, model.num_sec, model.n_run);

    case {'optimize','predict'}

        % free parameters
        aA                           = freeParam(1);
        bA                           = freeParam(2);
        sigV                         = freeParam(3);
        sigA                         = freeParam(4);
        sigP                         = freeParam(5);
        sigC                         = freeParam(6);
        pCommon                      = freeParam(7);

        sigMotor = model.sigMotor;
        muP = model.muP;

        if strcmp(model.mode, 'optimize')

            % constants
            CI.J_A            = sigA^2;
            CI.J_V            = sigV^2;
            CI.J_P            = sigP^2;
            CI.constC1        = CI.J_A*CI.J_V + CI.J_A*CI.J_P + CI.J_V*CI.J_P;
            CI.constC1_shat   = 1/CI.J_A  + 1/CI.J_V + 1/CI.J_P;
            CI.constC2_1      = CI.J_A + CI.J_P;
            CI.constC2_1_shat = 1/CI.J_A + 1/CI.J_P;
            CI.constC2_2      = CI.J_V + CI.J_P;
            CI.constC2_2_shat = 1/CI.J_V + 1/CI.J_P;

            % data: auditory locations (4) x visual locations (4) x postcues(2)
            % x rep
            [nLL_bimodal, ~] = calculateNLL_bimodal(...
                aA, bA, sigA, sigV, sigC, pCommon, ...
                CI, muP, sigMotor, data.bi_loc, data.bi_conf, model);

            out = sum(nLL_bimodal(:));

        elseif strcmp(model.mode, 'predict')

            out = sim_MA_gauss(aA, bA,  sigV, sigA, sigP, sigC, pCommon, model);

        end

end

end


function [nLL_bimodal, R] = calculateNLL_bimodal(...
    aA, bA, sigA, sigV, sigC, pCommon, ...
    CI, mu_P, sigMotor, data_resp, data_conf, model)

nLL_bimodal = 0;
sA_prime   = model.sA.*aA + bA; %the mean of biased auditory measurements
sV_prime  = model.sV;

n_sA = length(sA_prime);
n_sV = length(sV_prime);

for p = 1:length(sA_prime)   %for each AV pair with s_A' = s_A_prime(p)

    x1_grid    = linspace(sA_prime(p) - model.num_SD*sigA,...
        sA_prime(p) + model.num_SD*sigA,...
        model.numBins_A);
    mDist_AV_A = norm_dst(x1_grid, sA_prime(p), sigA, 0);


    for q = 1:length(sV_prime) %for each AV pair with s_V' = s_V_prime(q)

        x2_grid    = linspace(sV_prime(q) - model.num_SD*sigV,...
            sV_prime(q) + model.num_SD*sigV,...
            model.numBins_V);
        mDist_AV_V = norm_dst(x2_grid, sV_prime(q), sigV, 0);

        %use outer product to get the bi-variate distribution
        p_mAmV_given_sAsV = mDist_AV_V'.*mDist_AV_A;
        p_mAmV_given_sAsV = p_mAmV_given_sAsV./(sum(p_mAmV_given_sAsV(:)));

        %GRID_X1: a matrix where each row is a copy of x1_grid
        %GRID_X2: matrix where each column is a copy of x2_grid
        [GRID_X1, GRID_X2]= meshgrid(x1_grid, x2_grid);

        %calculate the likelihood and the posterior of a common cause and
        %separate causes L(C=1|x1=mA, x2=mV) and L(C=2|x1=mA, x2=mV)
        %size of Post_C1 (#rows x #cols): length(x2_grid) x length(x1_grid)
        [Post_C1, Post_C2, ~, ~]= calculatePostC1C2(GRID_X1, GRID_X2, CI, mu_P, pCommon);

        %calculate the shat_{C=1}, shat_{C=2} and MAP estimates assuming
        %model averaging
        shat_C1           = (GRID_X1./CI.J_A + GRID_X2./CI.J_V + ...
            mu_P./CI.J_P)./CI.constC1_shat;
        [shat_C2, MAP_MA] = deal(NaN(length(model.modality), model.numBins_V,...
            model.numBins_A));
        shat_C2(1,:,:) = (GRID_X1./CI.J_A + mu_P/CI.J_P)./CI.constC2_1_shat;
        shat_C2(2,:,:) = (GRID_X2./CI.J_V + mu_P/CI.J_P)./CI.constC2_2_shat;
        MAP_MA(1,:,:)  = shat_C1.*Post_C1 + squeeze(shat_C2(1,:,:)).*Post_C2;
        MAP_MA(2,:,:)  = shat_C1.*Post_C1 + squeeze(shat_C2(2,:,:)).*Post_C2;

        %--------------------- Confidence radius --------------------------

        % simulate posterior pdf for each trial using center coordinate
        post = zeros([2, model.numBins_V, model.numBins_A, numel(model.center_axis)]);
        for xx = 1:numel(model.center_axis)
            post(1,:,:,xx) = normpdf(model.center_axis(xx), squeeze(MAP_MA(1,:,:)), ...
                Post_C1/CI.constC1_shat + Post_C2/CI.constC2_1_shat + Post_C1*Post_C2*(shat_C1 - squeeze(shat_C2(1,:,:))).^2);
            post(2,:,:,xx) = normpdf(model.center_axis(xx), squeeze(MAP_MA(2,:,:)), ...
                Post_C1/CI.constC1_shat + Post_C2/CI.constC2_2_shat + Post_C1*Post_C2*(shat_C1 - squeeze(shat_C2(2,:,:))).^2);
        end

        % optimal radius given posterior and estimate
        post_2d = reshape(post, [prod([2, model.numBins_V, model.numBins_A]), numel(model.center_axis)]);
        shat_1d = reshape(MAP_MA, [prod([2, model.numBins_V, model.numBins_A]), 1]);
        [opt_radius, ~, ~]= eGain_MAP(post_2d, shat_1d, model.maxScore, model.minScore, model.elbow, model.center_axis);
        opt_radius = reshape(opt_radius, [2, model.numBins_V, model.numBins_A]);

        %----------------------- Compute likelihood -----------------------
        % For each same sA, sV combination, the data are organized by
        % response modality (A, V), repeititon of trials.
        % First half trials are A loc resps, second half trials are V loc
        % resps, from the same sA, sV combination
        locResp_A_V = [squeeze(data_resp(p,q,1,:))'; squeeze(data_resp(p,q,2,:))'];
        confResp_A_V = [squeeze(data_conf(p,q,1,:))'; squeeze(data_conf(p,q,2,:))'];
        [num_modality, num_rep] = size(locResp_A_V);

        for mm = 1:num_modality
            for kk = 1:num_rep
                % localization probability
                p_r_given_MAP = norm_dst(locResp_A_V(mm, kk), squeeze(MAP_MA(mm,:,:)),...
                    sigMotor,realmin);
                p_conf_given_m = norm_dst(confResp_A_V(mm, kk), squeeze(opt_radius(mm,:,:)),...
                    sigC, realmin);
                nLL_bimodal = nLL_bimodal - log(sum(sum(p_r_given_MAP.*...
                    p_conf_given_m.*p_mAmV_given_sAsV)));
            end
        end

        % %         elseif strcmp(model.strategy_loc,'MS') %Model selection
        %             for mm = 1:num_modality
        %                 for kk = 1:num_rep
        %                     % localization probability
        %                     if Post_C1(mm,kk) > 0.5
        %                         p_r_given_MAP = norm_dst(locResp_A_V(mm, kk),shat_C1,sigMotor,1e-20);
        %                     else
        %                         p_r_given_MAP = norm_dst(locResp_A_V(mm, kk),squeeze(shat_C2(mm,:,:)),sigMotor,1e-20);
        %                     end
        %                     % confidence probability
        %                     if confResp_A_V(mm, kk) == 1
        %                         p_conf_given_m = squeeze(p1_conf(mm,:,:));
        %                     elseif confResp_A_V(mm, kk) == 2
        %                         p_conf_given_m = squeeze(p2_conf(mm,:,:));
        %                     elseif confResp_A_V(mm, kk) == 3
        %                         p_conf_given_m = squeeze(p3_conf(mm,:,:));
        %                     elseif confResp_A_V(mm, kk) == 4
        %                         p_conf_given_m = squeeze(p4_conf(mm,:,:));
        %                     end
        %                     nLL_bimodal = nLL_bimodal - log(sum(sum(p_r_given_MAP.*...
        %                         p_conf_given_m.*p_mAmV_given_sAsV)));
        %                 end
        %             end

        %----------------------- Save if requested -----------------------
        R = [];
        if model.saveR

            %save the joint likelihood
            R.p_mAmV_given_sAsV(p,q,:,:) = p_mAmV_given_sAsV;

            %save predictions on location estimates
            R.loc(p,q,:,:,:) = MAP_MA;

            %             elseif strcmp(model.strategy_loc,'MS')
            %                 for mm = 1:length(model.modality)
            %                     R.loc(p,q,mm,:,:) = shat_C1;
            %                     R.loc(p,q,mm,Post_C1 < 0.5) = squeeze(shat_C2(mm,Post_C1 < 0.5));
            %                 end
            %             end

            %save predictions on confidence radius
            R.conf(p,q,:,:,:) = opt_radius;

        end

    end
end
end

function [Post_C1, Post_C2, L_C1, L_C2] = calculatePostC1C2(X1, X2, CI, mu_P, pC1)
%likelihood of a common cause and seperate causes
L_C1     = 1/(2*pi*sqrt(CI.constC1))*exp(-0.5*((X1 - X2).^2.*CI.J_P +...
    (X1 - mu_P).^2*CI.J_V + (X2 - mu_P).^2.*CI.J_A)./CI.constC1);
L_C2     = 1/(2*pi*sqrt(CI.constC2_1*CI.constC2_2))*exp(-0.5*...
    ((X1 - mu_P).^2./CI.constC2_1+(X2 - mu_P).^2./CI.constC2_2));
normTerm = L_C1.*pC1 + L_C2.*(1-pC1);
%posterior of a common cause
Post_C1  = L_C1.*pC1./normTerm;
Post_C2  = 1 - Post_C1;
end

function p = norm_dst(x,mu,sigma,t)
p = 1/sqrt(2*pi*sigma^2).*exp(-(x-mu).^2./(2*sigma^2)) + t;
end
