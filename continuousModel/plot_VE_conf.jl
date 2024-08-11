########## STANDARD LIBRARY PACKAGES ##########
# Pkg.add("ProgressBars", "ProgressMeter", "BenchmarkTools", "FileIO", "JLD2", "NPZ", "CSV", "DataFrames", "MAT", "JSON", "Plots", "LaTeXStrings", "StatsPlots")
# Pkg.add("CairoMakie")
using MAT
using Pkg
using Plots
using StatsPlots

# Change the working directory to the script's directory
script_dir = @__DIR__
cd(script_dir)

# Load the .mat file
model_names = ["Optimal", "MAP_optimal", "MAP_suboptimal", "Heuristic"]
sim_d = 2;
curr_model = model_names[sim_d];
filename = "sim_$curr_model.mat"
data = matread(filename)

# Access variables from the .mat file
for (key, value) in data
    # Use `@eval` to dynamically create variables
    @eval $(Symbol(key)) = $value
end
raw_diff = vec(raw_diff)
diffs = vec(diff)
n_cue = length(cue_label)
n_rel = length(rel_label)

# Plot VE
colors = [:brown4, :pink];  # Colors for rel=1 and rel=2
plt = plot(layout = (2, 2), size = (900, 700),
    titlefont = font("Helvetica", 14),
    guidefont = font("Helvetica", 12),
    tickfont = font("Helvetica", 10),
    legendfont = font("Helvetica", 12))
direction = [1, -1];

for cue in 1:n_cue
    for rel in 1:n_rel
        i_ve = ve_by_raw_diff[:, cue, rel]
        
        plot!(plt[cue], raw_diff, i_ve,
              color = colors[rel],
              linewidth = 2,
              xticks = round.(raw_diff), 
              tick_direction = :out,
              ylabel = "Shift of localization (cm)",
              xlabel = "Audiovisual discrepancy (V-A, cm)",
              title = cue_label[cue],
              label = (cue == 1 ? rel_label[rel] : false))  # Only label the first subplot
    end

    # y = 0
    plot!(plt[cue], raw_diff, zeros(size(raw_diff)),
          color = :black,
          linestyle = :solid,
          label = (cue == 1 ? "No effect" : false))  # Only label the first subplot
    
    # identity line
    plot!(plt[cue], raw_diff, direction[cue] .* raw_diff,
          linestyle = :dash,
          color = :black,
          label = (cue == 1 ? "full capture" : false))  # Only label the first subplot
end
plot!(plt[1], legendfontsize=8, legend=:bottomright) 

using StatsPlots

# Assuming plt and other variables (e.g., colors, rel_label, etc.) are already defined.

# Plot conf
ylim_min = minimum([minimum(conf_by_diff[i_diff][cue, rel, :]) for i_diff in 1:length(diffs), cue in 1:n_cue, rel in 1:n_rel])
ylim_max = maximum([maximum(conf_by_diff[i_diff][cue, rel, :]) for i_diff in 1:length(diffs), cue in 1:n_cue, rel in 1:n_rel])
sides = [:left, :right];  # Sides for rel=1 and rel=2
label_added = Dict(1 => false, 2 => false)
scaled_diffs = diffs / 6

for cue in 1:n_cue
    for rel in 1:n_rel
        i_conf = [vec(conf_by_diff[i_diff][cue, rel, :]) for i_diff in 1:length(diffs)]  # Gather all confidence data
        y_data = vcat(i_conf...)  # Flatten the list of vectors into a single vector
        x_data = vcat([fill(scaled_diffs[i_diff], length(conf_by_diff[i_diff][cue, rel, :])) for i_diff in 1:length(diffs)]...)

        # Use a ternary operator to set the label conditionally
        label = !label_added[rel] ? rel_label[rel] : raw""

        # Create the violin plot
        violin!(plt[cue+2], x_data, y_data,
            bandwidth = 1,
            show_median = true,
            side = sides[rel],  # Set the side of the violin plot
            linewidth = 0, 
            fillalpha = 0.7,
            color = colors[rel],  # Set the color based on rel
            xlims = (minimum(scaled_diffs) - 1, maximum(scaled_diffs) + 1),
            ylims = (0, ylim_max),
            tick_direction = :out,
            xticks = (scaled_diffs, Int.(round.(diffs))),
            ylabel = "Confidence radius (cm)",
            xlabel = "Audiovisual discrepancy (V-A, cm)",
            title = cue_label[cue],
            legendfontsize = 8,
            label = label)  # Conditionally add the label

        # Mark that the label has been added for this rel
        label_added[rel] = true
    end
end

display(plt)
