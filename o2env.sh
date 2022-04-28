
. /cvmfs/alice.cern.ch/etc/login.sh
eval $(alienv printenv VO_ALICE@QualityControl::v1.45.1-1)

# ALICE PyROOT is unhappy about the interference with cppyy from the standard
# SWAN environment, so we have to disable it:
unset CPPYY_BACKEND_LIBRARY


