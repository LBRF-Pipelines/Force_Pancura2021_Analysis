### MEP Preprocessing Pipeline Settings ###

# The sample rate (in Hz) of the recorded signal
srate <- 1000

# The onset (in ms) of the first TMS pulse.
pulse_on <- 5000

# The maximum possible TMS pulse onset (in ms) for the dataset
max_pulse_on <- 5000

# The window of time (in ms following each supra-threshold pulse) to search
# for an MEP. For example, for an SICI trial with a 120% RMT pulse at 1002 ms
# and a window of c(10, 90), the pipeline will look for the MEP between 1012 ms
# and 1092 ms, inclusive.
mep_window <- c(10, 50)

# The window of time (in ms) following a TMS pulse to look for the peak of the
# pulse artifact in the EMG signal
pulse_window <- 4

# The duration (in ms) of signal to keep for each frame. All samples after this
# timepoint will be discarded.
max_time <- 6000

# The frequency (in Hz) of the line noise to remove from the signal.
line_hz <- 60

# Whether or not to plot all processed MEPs as annotated PDFs for visual
# inspection.
plot_meps <- TRUE
