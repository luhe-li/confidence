function p_MAP = predict_prob_r(MAP_mat, p_mAmV_mat, sigma_r, model)
%This function generates predicted p(r|s_A, s_V) for all different
%combinations of s_A and s_V
%p_mAmV_mat: the joint probability of m_A and m_V given s_A and s_V
%            size = lenS x lenS x model.numBins_A x model.numBins_V
%MAP_mat   : the MAP estimates given an auditory and a visual stimulus
%            size = lenS x lenS x lenM x model.numBins_A x model.numBins_V
global lenS lenM
[p_r_given_shat,p_MAP] = deal(zeros(lenS, lenS, lenM, model.numBins_r));
for i = 1:lenS %for each auditory location
    for j = 1:lenS %for each visual location
        %get the joint likelihood and the matrix for MAP estimates given 
        %a selected visual and auditory location
        p_mAmV_temp = squeeze(p_mAmV_mat(i,j,:,:));
        MAP_temp = cell(1,lenM); 
        for k = 1:lenM; MAP_temp{k} = round(squeeze(MAP_mat(i,j,k,:,:)),1);end

        for l = 1:model.numBins_A 
            for m = 1:model.numBins_V
                for n = 1:lenM
                    %n = 1: get p(r_A|shat_A(l), shat_V(m))
                    %n = 2: get p(r_V|shat_A(l), shat_V(m))
                    p_r_given_shat(i,j,n,:) = norm_dst(model.bins_r,...
                        MAP_temp{n}(l,m),sigma_r,0);
                    p_MAP(i,j,n,:) = p_MAP(i,j,n,:) + ...
                        p_r_given_shat(i,j,n,:).*p_mAmV_temp(l,m);
                end
            end
        end
    end
end

function p = norm_dst(x,mu,sigma,t)
p = 1/sqrt(2*pi*sigma^2).*exp(-(x-mu).^2./(2*sigma^2)) + t;
