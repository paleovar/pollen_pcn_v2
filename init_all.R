# helper script to source all function definitions
#
# NOTE: when sourcing for the first time, the order matters 
# because some scripts depend on definitions in other scripts.
# Because of the handling of environments in R, re-definitions 
# which are actualized by Ctrl+Enter or Ctrl+Shift+S replace 
# existing ones. Thus, re-running `source("init_all.R")` with
# the order of initial sourcing as provided below is not 
# necessary after re-defining a function.

source('load_libraries.R')

source('constants.R')
source('util_dating.R')
source('util_file_processing.R')
source('util_signal_processing.R')
source('util_taxa_harmonization.R')
source('util_networks.R')
source('util_nullmodel.R')
source('util_plotting.R')
source('plot_maps.R')
source('plot_networks.R')
source('PCNdata.R')
