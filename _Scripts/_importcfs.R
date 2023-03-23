### Function for fast CFS import using reticulate ###

library(reticulate)
library(tibble)
library(purrr)


import_cfs <- function(f, chans, trialvars) {

    # Import pycfs module and raw CFS object
    pycfs <- import("pycfs")
    raw <- pycfs$read_cfs(f)

    # Extract trial metadata and build a dataframe
    n_frames <- length(raw$frames)
    meta_df <- tibble(frame = 1:n_frames)
    for (v in trialvars) {
        meta_df[[v]] <- pycfs$get_metadata(raw, v)
    }

    # Ensure consistent sample rates across data channels
    srates <- sapply(chans, function(ch) {
        raw$frames[[1]]$sample_rates[[ch]]
    })
    if (length(unique(srates)) > 1) {
        stop("Unable to parse multiple sample rates per frame.")
    }
    srate <- srates[1]
    
    # Extract signal data and build a dataframe
    signal_df <- map_df(1:n_frames, function(n) {
        n_samples <- length(raw$frames[[n]]$data[[chans[[1]]]])
        framedat <- tibble(
            frame = n,
            time = 1:n_samples * srate
        )
        for (ch in chans) {
            framedat[[ch]] <- raw$frames[[n]]$data[[ch]]
        }
        framedat
    })
    
    list(trial = meta_df, signal = signal_df)
}