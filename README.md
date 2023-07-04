# Grip Force and Motor Imagery - Pancura 2021

Analysis code for a study examining the difference in corticospinal excitability between imagined and overt movements performed at varying forces.

## Dependencies

To install the R packages required for running this analyis pipeline, run the following command at an R prompt:

```r
install.packages(
  c("tidyverse", "arrow", "pracma", "brms", "bayestestR" "emmeans")
)
```

## Running The Pipeline

First, to run the pipeline, you will need to place the contents of the `Data` folder from the project's OSF repository (found [here](https://osf.io/z8ct6/files/osfstorage)) in the pipeline's `_Data` folder.

Then, set the working directory in R to the `Force_Pancura2021_Analysis` folder and run the following source commands in sequence:

```r
# Imports all task & signal data for the study
source("./_Scripts/0_import.R")

# Preprocesses & summarizes force transducer data
source("./_Scripts/1_force.R")

# Preprocesses & summarizes MEP and background EMG data
source("./_Scripts/2_processing.R")

# Performs final data cleaning and runs mixed-effects models
source("./_Scripts/3_modelling.R")

# Summarizes demographics, exclusions, and model results
source("./_Scripts/4_reporting.R")

# Generates plots based on data
source("./_Scripts/5_visualization.R")

```

