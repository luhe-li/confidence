# install.packages("R.matlab")
# install.packages("patchwork")

library(R.matlab)
library(rstudioapi)
library(ggplot2)
library(patchwork)

# Set working directory to the location of the sourcing file
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
D <- readMat('sim_data.mat')

# Define variables
model_names <- c('Optimal', 'MAP optimal', 'MAP suboptimal', 'Heuristic')
cue_label <- c('Post-cue: A', 'Post-cue: V')
rel_label <- c('High visual reliability', 'Low visual reliability')

# Extract variables
sim_d <- D$sim.d
n_rep <- D$n.rep
raw_diff <- D$raw.diff
ve_by_raw_diff <- D$ve.by.raw.diff
model_name <- model_names[sim_d]
n_cue <- length(cue_label)
n_rel <- length(rel_label)

# Prepare the plot title
main_title <- sprintf("%s, rep: %i", model_name, n_rep)
# Create a list to store the plots
plots <- list()

for (cue in 1:n_cue) {
  p <- ggplot() +
    labs(title = cue_label[cue],
         x = "Audiovisual discrepancy (V-A, cm)",
         y = "Shift of localization (cm)") +
    theme_minimal() +
    coord_fixed() +
    # theme(plot.title = element_text(hjust = 0.5),
    #       panel.grid = element_blank(),
    #       plot.margin = margin(10, 10, 10, 10)) +
    xlim(min(raw_diff), max(raw_diff)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
  
  for (rel in 1:dim(ve_by_raw_diff)[3]) {
    i_ve <- ve_by_raw_diff[, cue, rel]
    p <- p + geom_line(aes(x = raw_diff, y = i_ve))
  }
  
  if (cue == 1) {
    p <- p + geom_line(aes(x = raw_diff, y = raw_diff), linetype = "dashed", color = "black")
  } else {
    p <- p + geom_line(aes(x = raw_diff, y = -raw_diff), linetype = "dashed", color = "black")
  }
  
  plots[[cue]] <- p
}

# Combine the plots into a 1x2 layout
combined_plot <- plots[[1]] + plots[[2]] + plot_layout(ncol = 2)

# Display the combined plot
print(combined_plot)