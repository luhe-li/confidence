function [pred_mean_conf, abs_diff] = analyze_pred(pred)

% organize confidence data: {diff} cue x reliability x rep
[conf_by_diff, abs_diff] = org_by_diffs(pred.bi_conf, pred.sA);

for cue = 1:2
    for rel = 1: 2
        for diff = 1:numel(abs_diff)
            i_conf = squeeze(conf_by_diff{diff}(cue, rel, :))';
            p = sum(i_conf)/(numel(i_conf));
            pred_mean_conf(diff, cue, rel) = p;
        end
    end
end



end