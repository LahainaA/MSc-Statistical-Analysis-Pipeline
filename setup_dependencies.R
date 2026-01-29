# SETUP: INSTALL DEPENDENCIES

required_packages <- c(
  "tidyverse",    # Includes ggplot2, dplyr, readr, tidyr, etc.
  
  # Statistics & Modeling
  "lme4",         # Linear Mixed Effects Models
  "lmerTest",     # P-values for LMEMs
  "afex",         # ANOVA convenience wrappers 
  "emmeans",      # Post-hoc contrasts and marginal means
  "minpack.lm",   # Non-linear least squares (Required for Landmark/Sine fit)
  
  # Diagnostics & Validation
  "DHARMa",       # Residual diagnostics for LMEMs
  
  # Reporting & Visualization
  "knitr",        # R Markdown report generation
  "rmarkdown",    # Compiling .Rmd files
  "kableExtra",   # Pretty tables (Demographics)
  "patchwork",    # Combining plots 
  "gridExtra"     # Combining plots 
)

# Check and Install
print("Checking for missing packages...")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if(length(missing_packages)) {
  print(paste("Installing:", paste(missing_packages, collapse = ", ")))
  install.packages(missing_packages)
} else {
  print("All packages are already installed.")
}

print("Setup Complete! You can now run the .Rmd files.")