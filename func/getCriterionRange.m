function [lb, ub] = getCriterionRange(sigV1, sigA, sigV2, sigP, model)

% This function calculates a reasonable range of criterion of variance
% given the noise of auditory, high and low visual stimulus

% Based on previous simulation, the optimal and suboptimal model
% can share the same criterion range

JA = sigA^2;
JV1 = sigV1^2;
JV2 = sigV2^2;
JP = sigP^2;

switch model

    case {'Optimal','Suboptimal'}
        % suboptimal and optimal
        ub = max([1/(1/JV1 + 1/JA + 1/JP), 1/(1/JV2 + 1/JA + 1/JP)]);
        lb = min([1/(1/JV1 + 1/JP), 1/(1/JV2 + 1/JP)]);

    case 'Heuristic'
        % heuristic
        ub = max([JA, JV1, JV2]);
        lb = min([JA, JV1, JV2]);

end