get_map_from_samples <- function(samples, med = TRUE, breaks = 10, plot = FALSE, ...) {
    ## get hist 
    density_x <- hist(samples, breaks = breaks, plot = plot, ...)
    ## get max density
    max_density <- max(density_x$density)
    ## decrease the density until 95% of sample are covered
    
    idx_max <- which.max(density_x$density)
    if (med) {
        median(samples[samples > density_x$breaks[idx_max] &
                         samples < density_x$breaks[idx_max + 1]])
    } else{
        mean(samples[samples > density_x$breaks[idx_max] &
                         samples < density_x$breaks[idx_max + 1]])
    }
}