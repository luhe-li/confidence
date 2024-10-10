function out = nll_uni(free_param, model, data)

if strcmp(model.mode, 'initiate')

    out.param_id = {'aA','bA','\sigma_{A}','\sigma_{V]','\sigma_C','\sigma_{P}','\mu_P'};
    out.num_param = length(out.param_id);

    % hard bounds, the range for lb, ub, larger than soft bounds
    param_h.aA = [-3, 3];
    param_h.bA = [-10, 10];
    param_h.sigA = [1, 20];
    param_h.sigV = [1e-2, 5];
    param_h.sigC = [1e-2, 10];
    param_h.sigP = [0.01, 10];
    param_h.muP = [-5, 5];

    % soft bounds, the range for plb, pub
    param_s.aA = [-1.5, 1.5];
    param_s.bA = [-5, 5];
    param_s.sigA = [5, 10];
    param_s.sigV = [1, 3];
    param_s.sigC = [1, 5];
    param_s.sigP = [1, 5];
    param_s.muP = [-1, 1];

    % reorganize parameter bounds to feed to bads
    fields = fieldnames(param_h);
    for k = 1:numel(fields)
        out.lb(:,k) = param_h.(fields{k})(1);
        out.ub(:,k) = param_h.(fields{k})(2);
        out.plb(:,k) = param_s.(fields{k})(1);
        out.pub(:,k) = param_s.(fields{k})(2);
    end
    model.param_s = param_s;
    model.param_h = param_h;

    % get grid initializations
    num_sections = model.n_run*2;
    out.init = getInit(out.lb, out.ub, num_sections, model.n_run);

else

    % assign free parameters
    aA = free_param(1);
    bA = free_param(2);
    sigma_A = free_param(3);
    sigma_V = free_param(4);
    sigma_C = free_param(5);
    sigma_P = free_param(6);
    mu_P = free_param(7);

    if strcmp(model.mode, 'optimize')

        %------------------------ localization-----------------------------
        s_A_prime_uni = repmat(model.uni_sA',[1, model.uni_nrep]).* aA + bA;
        s_V_prime_uni = repmat(model.uni_sV',[1, model.uni_nrep]) .* model.aV + model.bV;

        %calculate the mean of the response distributions
        c_A                = (1/sigma_A^2)/(1/sigma_A^2+1/sigma_P^2);
        c_V                = (1/sigma_V^2)/(1/sigma_V^2+1/sigma_P^2);
        f_A                = (mu_P/sigma_P^2)/(1/sigma_A^2+1/sigma_P^2);
        f_V                = (mu_P/sigma_P^2)/(1/sigma_V^2+1/sigma_P^2);

        %the means of estimate distributions
        mu_shat_A_uni      = reshape(c_A.*s_A_prime_uni + f_A, [1, numel(model.uni_sA)*model.uni_nrep]);
        mu_shat_V_uni      = reshape(c_V.*s_V_prime_uni + f_V, [1, numel(model.uni_sV)*model.uni_nrep]);

        %the variances of estimate distributions
        sigma_shat_A       = c_A*sigma_A;
        sigma_shat_V       = c_V*sigma_V;

        % variance of the localization response distributions
        sigma_a_loc_resp = sqrt(sigma_shat_A^2 + model.sigma_motor^2);
        sigma_v_loc_resp = sqrt(sigma_shat_V^2 + model.sigma_motor^2);

        % reshape data into 2 rows (1st row: A localization responses,
        % location and nrep flattened; 2end row: V localization responses)
        uni_loc_2d = reshape(permute(data.uni_loc, [2, 1, 3]), model.modality, []);
        uni_conf_2d = reshape(permute(data.uni_conf, [2, 1, 3]), model.modality, []);

        nLL_loc_unimodal    = calculateNLL_unimodal([mu_shat_A_uni; mu_shat_V_uni], ...
            [sigma_a_loc_resp; sigma_v_loc_resp], uni_loc_2d);

        %------------------------ confidence -----------------------------
        % simulate optimal confidence radius, which is independent of
        % measurement
        sim = sim_uni(aA, bA, sigma_A, sigma_V, sigma_P, mu_P, sigma_C, model);
        opt_radius = reshape(permute(sim.opt_radius, [2, 1, 3]), model.modality, []);

        nLL_conf_unimodal    = calculateNLL_unimodal(opt_radius, ...
            [sigma_C; sigma_C], uni_conf_2d);

        out = nLL_loc_unimodal+nLL_conf_unimodal;

    elseif strcmp(model.mode, 'predict')

        out = sim_uni(aA, bA, sigma_A, sigma_V, sigma_P, mu_P, sigma_C, model);

    end
end

    function nLL_unimodal = calculateNLL_unimodal(mu, sig, x)

        % --------------------- Localization ------------------------------
        %we assume that response distributions are Gaussian, centered at mu
        %with variance equal to (sigma_shat^2 + sigma_motor^2)
        %mu, sig, x all consist of 2 rows (1st row: A; 2nd row: V)

        LL = arrayfun(@(idx) length(mu(idx,:))*(-0.5*log(2*pi*sig(idx)^2))-...
            sum((x(idx,:) - mu(idx,:)).^2)./(2*sig(idx)^2), 1:2);
        nLL_unimodal = sum(-LL(:));

    end


end