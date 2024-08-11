########## STANDARD LIBRARY PACKAGES ##########
# Pkg.add("ProgressBars", "ProgressMeter", "BenchmarkTools", "FileIO", "JLD2", "NPZ", "CSV", "DataFrames", "MAT", "JSON", "Plots", "LaTeXStrings", "StatsPlots")
using MAT
using Pkg
using Plots
using StatsPlots
using GLMakie
using CairoMakie

# Change the working directory to the script's directory
script_dir = @__DIR__
cd(script_dir)

# Load the .mat file
file = matread("sim_data.mat")

# Access variables from the .mat file
for (key, value) in file
    # Use `@eval` to dynamically create variables
    @eval $(Symbol(key)) = $value
end
raw_diff = vec(raw_diff)
diffs = vec(diff)
n_cue = length(cue_label)
n_rel = length(rel_label)

# Plot VE
colors = [:brown4, :pink];  # Colors for rel=1 and rel=2
plt = plot(layout = (2, 2), size = (800, 600))

for cue in 1:n_cue
    for rel in 1:n_rel
        i_ve = ve_by_raw_diff[:, cue, rel];
        plot!(plt[cue], raw_diff, i_ve,
        color = colors[rel],
        linewidth = 2,
        xticks = round.(raw_diff), 
        tick_direction = :out,
        ylabel = "Shift of localization (cm)",
        xlabel = "Audiovisual discrepancy (V-A, cm)",
        title = cue_label[cue])
    end
    
    if cue == 1
        plot!(plt[cue], raw_diff, raw_diff, linestyle = :dash, color = :black)
    else
        plot!(plt[cue], raw_diff, -raw_diff, linestyle = :dash, color = :black)
    end
end

# Plot conf

ylim_min = minimum([minimum(conf_by_diff[i_diff][cue, rel, :]) for i_diff in 1:length(diff), cue in 1:n_cue, rel in 1:n_rel])
ylim_max = maximum([maximum(conf_by_diff[i_diff][cue, rel, :]) for i_diff in 1:length(diffs), cue in 1:n_cue, rel in 1:n_rel])
sides = [:left, :right];  # Sides for rel=1 and rel=2
label_added = Dict(1 => false, 2 => false);

for cue in 1:n_cue
    for rel in 1:n_rel
        i_conf = [vec(conf_by_diff[i_diff][cue, rel, :]) for i_diff in 1:length(diffs)];  # Gather all confidence data

        # Convert to a flat vector and generate corresponding x positions with jitter
        y_data = vcat(i_conf...);  # Flatten the list of vectors into a single vector

        # Use a ternary operator to set the label conditionally
        label = !label_added[rel] ? rel_label[rel] : raw""

        # Create the violin plot
        violin!(plt[cue+2], diffs, y_data,
            bandwidth = 0.5,
            # scale = :count,
            side = sides[rel],  # Set the side of the violin plot
            linewidth = 0, 
            fillalpha = 0.5,
            color = colors[rel],  # Set the color based on rel
            xlims = (minimum(diffs),maximum(diffs)),  # Adjusted x limits to fit the new xticks
            ylims = (ylim_min, ylim_max),
            tick_direction = :out,
            xticks = (round.(diffs)),  # Set xticks positions and corresponding labels from diffs
            ylabel = "Confidence radius (cm)",
            xlabel = "Audiovisual discrepancy (V-A, cm)",
            title = cue_label[cue],
            label = label)  # Conditionally add the label

        # Mark that the label has been added for this rel
        label_added[rel] = true
    end
end

display(plt)
###

# ylim_min = minimum([minimum(conf_by_diff[i_diff][cue, rel, :]) for i_diff in 1:length(diffs), cue in 1:n_cue, rel in 1:n_rel])
# ylim_max = maximum([maximum(conf_by_diff[i_diff][cue, rel, :]) for i_diff in 1:length(diffs), cue in 1:n_cue, rel in 1:n_rel])
# sides = [:left, :right];  # Sides for rel=1 and rel=2
# label_added = Dict(1 => false, 2 => false);

# for cue in 1:n_cue
#     for rel in 1:n_rel
#         all_conf = [vec(conf_by_diff[i_diff][cue, rel, :]) for i_diff in 1:length(diffs)]  # Gather all confidence data

#         # Convert to a flat vector and generate corresponding x positions with jitter
#         y_data = vcat(all_conf...);  # Flatten the list of vectors into a single vector

#         # Use a ternary operator to set the label conditionally
#         label = !label_added[rel] ? rel_label[rel] : ""

#         violin!(plt[cue+2], diff, y_data,
#             bandwidth = 1,
#             side = sides[rel],  # Set the side of the violin plot
#             linewidth = 0, fillalpha = 0.5,
#             color = colors[rel],  # Set the color based on rel
#             xlim = (minimum(diff) - 5, maximum(diff) + 5),
#             ylims = (ylim_min, ylim_max),
#             tick_direction = :out,
#             xticks = round.(diff),
#             ylabel = "Confidence radius (cm)",
#             xlabel = "Audiovisual discrepancy (V-A, cm)",
#             title = cue_label[cue],
#             label = label)  # Conditionally add the label

#         # Mark that the label has been added for this rel
#         label_added[rel] = true
#     end
# end

# # Display the plot
# display(plt)
