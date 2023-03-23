##########################################
### Utility Functions for MEP Analysis ###
##########################################

library(dplyr)
library(tibble)
library(pracma)


### Utility Functions ###

# Prints output to console if the variable "verbose" exists and is TRUE

optlog <- function(msg) {
  if ("verbose" %in% ls(envir = .GlobalEnv) && verbose == TRUE) {
    cat(msg)
  }
}


# NOTE: this is copied almost directly from the 'pracma' package. The actual
# function won't work because it mistakenly errors on complex numbers, so
# this modified version fixes that.

pinv <- function (A, tol = .Machine$double.eps ^ (2 / 3)) {

  stopifnot(
    is.numeric(A) | is.complex(A), length(dim(A)) == 2, is.matrix(A)
  )
  s <- svd(A)

  p <- ( s$d > max(tol * s$d[1], 0) )
  if (all(p)) {
      mp <- s$v %*% (1 / s$d * t(s$u))
  } else if (any(p)) {
      mp <- s$v[, p, drop = FALSE] %*% (1 / s$d[p] * t(s$u[, p, drop = FALSE]))
  } else {
      mp <- matrix(0, nrow = ncol(A), ncol = nrow(A))
  }

  mp
}


### Signal Processing Functions ###

notch_filter <- function(signal, freq, srate) {

  # Create Pei-Tseng notch filter
  nyquist_hz <- srate / 2
  notch <- gsignal::butter(4, c(freq - 2, freq + 2) / nyquist_hz, "stop")

  # Apply notch filter and return filtered data
  gsignal::filtfilt(notch, signal)
}


# A line noise filter that works really well, but causes post-MEP distortions
# for higher-amplitude MEPs (a problem for LICI trials). Included here for
# comparitive purposes.

pei_tseng_filter <- function(signal, freq, srate) {

  # Create Pei-Tseng notch filter
  nyquist_hz <- srate / 2
  notch <- gsignal::pei_tseng_notch(freq / nyquist_hz, 2 / nyquist_hz)

  # Apply notch filter and return filtered data
  gsignal::filter(notch, signal)
}


# A line noise filter based on discrete Fourier transformation (DFT).
# Based on the 'zero' method of FieldTrip's ft_preproc_dftfilter, except it
# only supports a single channel at a time.

dft_filt <- function(signal, line_hz, srate, window = NA) {

  # Get the maximum number of samples divisible by the line noise frequency
  nsamples <- length(signal)
  if (is.na(window)) {
    n <- as.integer(nsamples / line_hz) * line_hz
  } else {
    n <- as.integer(window / line_hz) * line_hz
  }

  # Temporary subtract the mean from the signal
  sigmean <- mean(signal)
  signal <- signal - sigmean

  # Generate a complex sine wave at the line noise frequency
  time <- (1:nsamples - 1) / srate
  k <- 2 * pi * line_hz * time
  noise_tmp <- exp(1i * k)

  # Estimate the amplitude of the line noise from the signal
  ampl <- 2 * t(signal[1:n]) %*% pinv(t(noise_tmp[1:n]))
  noise_est <- noise_tmp * as.vector(ampl)

  # Subtract the estimated line noise from the signal and re-add the mean
  filtered <- Re(signal - noise_est) + sigmean

  filtered
}


dft_filt_sine <- function(signal, line_hz, srate, window = NA) {
  
  # Get the maximum number of samples divisible by the line noise frequency
  nsamples <- length(signal)
  if (is.na(window)) {
    n <- as.integer(nsamples / line_hz) * line_hz
  } else {
    n <- as.integer(window / line_hz) * line_hz
  }
  
  # Temporary subtract the mean from the signal
  sigmean <- mean(signal)
  signal <- signal - sigmean
  
  # Generate a complex sine wave at the line noise frequency
  time <- (1:nsamples - 1) / srate
  k <- 2 * pi * line_hz * time
  noise_tmp <- exp(1i * k)
  
  # Estimate the amplitude of the line noise from the signal
  ampl <- 2 * t(signal[1:n]) %*% pinv(t(noise_tmp[1:n]))
  noise_est <- noise_tmp * as.vector(ampl)
  
  # Subtract the estimated line noise from the signal and re-add the mean
  sine <- Re(noise_est) + sigmean
  
  sine
}


# A function for upsampling signals using spline interpolation.

spline_upsample <- function(signal, time, new_srate) {

  # Calculate the current and new periods for the signal
  current_period <- min(time - lag(time), na.rm = TRUE)
  new_period <- 1 / new_srate

  # If the new sample rate is higher, do spline interpolation
  if (new_period < current_period) {
    f <- splinefun(time, signal, method = "monoH.FC")
    time <- seq(min(time), max(time), by = new_period)
    signal <- f(time)
  }

  tibble(time = time, signal = signal)
}



### Peak Detection Functions ###

detect_peaks <- function(signal, time, min_srate) {

  # NOTE: ignoring secondary peaks?

  # Look for MEP min/max values within the window
  mep_max <- max(signal)
  mep_min <- min(signal)

  # Get timestamps of MEP max/min values
  mep_max_t <- min(time[signal == mep_max])
  mep_min_t <- min(time[signal == mep_min])

  # Return list of values
  tibble(
    vmax = mep_max, vmin = mep_min,
    tmax = mep_max_t, tmin = mep_min_t
  )

}


detect_onset <- function(signal, time, baseline_var) {

  near_baseline <- signal < baseline_var * 2
  before_spike <- near_baseline & !lead(near_baseline)

  if (any(before_spike, na.rm = TRUE)) {
    onset <- max(which(before_spike))
  } else {
    onset <- length(time) - 1  # default to sample before peak
  }

  time[onset]
}



### Visualization Functions ###

plot_psd <- function(signal, srate, db = TRUE) {
  psd_dat <- spec.pgram(signal, pad = 1, taper = 0.2, plot = FALSE)
  freqs_hz <- psd_dat$freq * srate
  pwr <- psd_dat$spec
  if (db) {
    pwr <- log(pwr)
  }
  ggplot() +
    geom_line(aes(x = freqs_hz, y = pwr)) +
    xlab("Frequency (Hz)") +
    ylab("Power (dB)")
}

plot_mep <- function(amp, time, window = c(-0.5, 0.2)) {
  ggplot() +
    geom_line(aes(x = time, y = amp)) +
    geom_vline(x_intercept = peaks)
}