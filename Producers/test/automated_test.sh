# bash script

#####################################################################
# Automates scram, batch_condor.py and cmsRun of a test sample      #
# This is meant to be run interactively, e.g. `. automated_test.sh` #
#                                                                   #
# Author: A. Mecca (amecca@cern.ch)                                 # 
#####################################################################

set -o pipefail

analyzer=../python/analyzer_VVjj.py  # ../python/analyzer_VVjj_dump.py  # 
specific_file=""  # es. file:/eos/home-a/amecca/samples/.../somefile.root
csv_file=samples_2017UL_MC.csv  # ../python/samples_PKU_2016UL_MC.csv
sample_ID=ZGToLLG  # ZGToLLG-v17
events=100


# First, recompile CMSSW and remove any results from last execution
cd $CMSSW_BASE/src ; printf "scram... " ; scram b -j 12 1>/dev/null && echo "Done" || { cd - ; echo ; return 1 2>/dev/null || exit 1 ; }
cd $CMSSW_BASE/src/VVXAnalysis/Producers/test || { echo "No VVXAnalysis in CMSSW/src" 1>&2 ; return 1 2>/dev/null || exit 1 ; }
[ -d auto_test ] && { rm -r auto_test || { echo "Failed to rm auto_test" 1>&2 ; return 1 2>/dev/null || exit 1 ; } ; }

# Create the config using the specified template ("analyzer") for all the uncommented lines in csv_file
batch_Condor.py ${csv_file} -i $analyzer -o auto_test | grep -Pv "N:\s+\d*[1-9]\s" || { return 1 2>/dev/null || exit 1 ; }

# Check which samples were created
samples=$(ls auto_test/ | grep -oP ".+(?=run_template_cfg.py)")
[ -z "$samples" ] && { echo "ERROR: no samples were prepared" 1>&2 ; return 1 2>/dev/null || exit ; }
if ! echo $samples | grep -q $sample_ID ; then
    printf 'WARN: "%s" was not in the csv. ' $sample_ID
    sample_ID=$(echo "$samples" | head -n1)
    printf 'Using "%s"\n' $sample_ID
fi

# Take the first chunk, which should always exist
chunk=auto_test/${sample_ID}_Chunk0
cd $chunk && echo "INFO: running in $chunk" || { echo "cd into $chunk failed." 1>&2 ; return 1 2>/dev/null || exit 1 ; }

# Set the maximum number of events
printf "process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(%d) )\n" $events >> run_cfg.py

# Optional: insert a custom filename
[ -n "$specific_file" ] && -i -r "s/(fileNames = cms\.untracked\.vstring).+$/\1('${specific_file}'),/" run_cfg.py 1>/dev/null

# Run in background, redirect stdout and stderr to a file and read from it (so that the log is permanent and you can ^C without killing cmsRun)
cmsRun run_cfg.py &> cmsrun.log &
less +F cmsrun.log
cd ../..
