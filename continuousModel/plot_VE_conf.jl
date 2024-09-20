########## STANDARD LIBRARY PACKAGES ##########
# Pkg.add("ProgressBars", "ProgressMeter", "BenchmarkTools", "FileIO", "JLD2", "NPZ", "CSV", "DataFrames", "MAT", "JSON", "Plots", "LaTeXStrings", "StatsPlots")
# Pkg.add("CairoMakie")
using MAT
using Pkg
using Plots
using Statistics  # Add this line to import the Statistics package

# Change the working directory to the script's directory
script_dir = @__DIR__
cd(script_dir)

# Load the .mat file
model_names = ["Optimal model average", "Lazy model average","Model selection","Probability matching"];
folders = ["optimal","MA","MS","PM"];
sim_d = 2;
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
elbow_value = round(fixP["elbow"], digits = 1)

# Define custom colors
clt = [5 113 176; 202 0 32] ./ 255
colors = [RGB(clt[1, 1], clt[1, 2], clt[1, 3]), RGB(clt[2, 1], clt[2, 2], clt[2, 3])]  # Colors for cue=1 and cue=2

# Plot VE
l = @layout [a{0.01h}; grid(2,2)]
plt = plot(layout = l, size = (900, 700), dpi=300,
    titlefont = font("Helvetica", 14),
    guidefont = font("Helvetica", 12),
    tickfont = font("Helvetica", 10),
    legendfont = font("Helvetica", 8))
i_model_name = model_names[Int(sim_d)]
plot!(plt[1], title= "$i_model_name, elbow = $elbow_value cm", framestyle=nothing, showaxis=false, xticks=false, yticks=false)
direction = [1, -1];

for cue in 1:n_cue
    i_ve = ve_by_raw_diff[:, cue, 1]  # Only plot for rel=1
    
    plot!(plt[cue+1], raw_diff, i_ve,
          color = colors[cue],
          linewidth = 2,
          xticks = round.(raw_diff), 
          tick_direction = :out,
          ylims = [-15, 15],
          ylabel = "Shift of localization (cm)",
          xlabel = "Audiovisual discrepancy (V-A, cm)",
          title = cue_label[cue],
          label = (cue == 1 ? rel_label[1] : false),
          legend = :bottomright)  # Only label the first subplot

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
ylim_min = minimum([minimum(conf_by_diff[i_diff][cue, 1, :]) for i_diff in 1:length(diffs), cue in 1:n_cue])
ylim_max = maximum([maximum(conf_by_diff[i_diff][cue, 1, :]) for i_diff in 1:length(diffs), cue in 1:n_cue])
scaled_diffs = diffs / 10

for cue in 1:n_cue
    i_conf = [vec(conf_by_diff[i_diff][cue, 1, :]) for i_diff in 1:length(diffs)]  # Gather all confidence data
    y_data = vcat(i_conf...)  # Flatten the list of vectors into a single vector
    x_data = vcat([fill(scaled_diffs[i_diff], length(conf_by_diff[i_diff][cue, 1, :])) for i_diff in 1:length(diffs)]...)

    # Create the violin plot with overlayed box plot
    violin!(plt[cue+3], x_data, y_data,
        bandwidth = 1,
        show_median = true,
        side = :both,  # Full violin plot
        linewidth = 0, 
        fillalpha = 1,
        color = colors[cue],  # Set the color based on cue
        xlims = (minimum(scaled_diffs) - 1, maximum(scaled_diffs) + 1),
        ylims = (0, 25),
        tick_direction = :out,
        xticks = (scaled_diffs, Int.(round.(diffs))),
        ylabel = "Confidence radius (cm)",
        xlabel = "Audiovisual discrepancy (V-A, cm)",
        title = cue_label[cue],
        label = false)
    
    # Overlay median line
    for i_diff in 1:length(diffs)
        median_value = median(conf_by_diff[i_diff][cue, 1, :])
        plot!(plt[cue+3], [scaled_diffs[i_diff] - 0.4, scaled_diffs[i_diff] + 0.4], [median_value, median_value],
              color = :black,
              linewidth = 3,
              label = false)
    end
end
display(plt)

a = RGBA(0.75, 0.75, 0.75, 0.5)
# save figure
# Step 1: Get the name of the current script
script_name = basename(@__FILE__)  # Get the script name with extension
script_name_without_extension = splitext(script_name)[1]  # Remove the extension

# Step 2: Create a directory with the script name if it doesn't exist
output_dir = joinpath(pwd(), script_name_without_extension)  # Full path in the current directory
if !isdir(output_dir)
    mkpath(output_dir)
end
savefig(joinpath(output_dir, "sim_$curr_model elbow $elbow_value.png"))
