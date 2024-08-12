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
model_names = ["Optimal", "Model average","Model selection","Heuristic"];
folders = ["optimal","MA","MS","Heuristic"];
sim_d = 1;
curr_model = folders[sim_d];
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
l = @layout [a{0.01h}; grid(2,2)]
plt = plot(layout = l, size = (900, 700), dpi=300,
    titlefont = font("Helvetica", 14),
    guidefont = font("Helvetica", 12),
    tickfont = font("Helvetica", 10),
    legendfont = font("Helvetica", 8))
plot!(plt[1], title=model_names[Int(sim_d)],framestyle=nothing,showaxis=false,xticks=false,yticks=false)
direction = [1, -1];

for cue in 1:n_cue
    for rel in 1:n_rel
        i_ve = ve_by_raw_diff[:, cue, rel]
        
        plot!(plt[cue+1], raw_diff, i_ve,
              color = colors[rel],
              linewidth = 2,
              xticks = round.(raw_diff), 
              tick_direction = :out,
              ylabel = "Shift of localization (cm)",
              xlabel = "Audiovisual discrepancy (V-A, cm)",
              title = cue_label[cue],
              label = (cue == 1 ? rel_label[rel] : false),
              legend = :bottomright)  # Only label the first subplot
    end

    # y = 0
    plot!(plt[cue+1], raw_diff, zeros(size(raw_diff)),
          color = :black,
          linestyle = :solid,
          label = (cue == 1 ? "No effect" : false))  # Only label the first subplot
    
    # identity line
    plot!(plt[cue+1], raw_diff, direction[cue] .* raw_diff,
          linestyle = :dash,
          color = :black,
          label = (cue == 1 ? "Full capture" : false))  # Only label the first subplot
end

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
        violin!(plt[cue+3], x_data, y_data,
            bandwidth = 1,
            show_median = true,
            side = sides[rel],  # Set the side of the violin plot
            linewidth = 0, 
            fillalpha = 1,
            color = colors[rel],  # Set the color based on rel
            xlims = (minimum(scaled_diffs) - 1, maximum(scaled_diffs) + 1),
            ylims = (0, 21),
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

# save figure
# Step 1: Get the name of the current script
script_name = basename(@__FILE__)  # Get the script name with extension
script_name_without_extension = splitext(script_name)[1]  # Remove the extension

# Step 2: Create a directory with the script name if it doesn't exist
output_dir = joinpath(pwd(), script_name_without_extension)  # Full path in the current directory
if !isdir(output_dir)
    mkpath(output_dir)
end
savefig(joinpath(output_dir, "sim_$curr_model.png"))