#!/bin/env python3

# Creates HTCondor jobs for an EventAnayzer

import os
import stat
import logging
from subprocess import check_call, check_output #, CalledProcessError
from argparse import ArgumentParser


condorsub_template = '''\
executable              = $(directory)/batchScript.sh
arguments               = {mainDir}/$(directory) $(ClusterId)$(ProcId)
output                  = log/$(ClusterId).$(ProcId).out
error                   = log/$(ClusterId).$(ProcId).err
log                     = log/$(ClusterId).$(ProcId).log
Initialdir              = $(directory)
request_memory          = 4000M
#Possible values: https://batchdocs.web.cern.ch/local/submit.html
+JobFlavour             = "{jobFlavour}"

x509userproxy           = {home}/x509up_u{uid}

#https://www-auth.cs.wisc.edu/lists/htcondor-users/2010-September/msg00009.shtml
periodic_remove         = JobStatus == 5

ShouldTransferFiles     = YES
'''

script_template = '''\
#!/bin/bash
set -euox pipefail

if [ -z ${{_CONDOR_SCRATCH_DIR+x}} ]; then
  #running locally
  runninglocally=true
  _CONDOR_SCRATCH_DIR=$(mktemp -d)
  SUBMIT_DIR=$(pwd)
else
  runninglocally=false
  SUBMIT_DIR=$1
fi

cd $SUBMIT_DIR
make

cd $_CONDOR_SCRATCH_DIR

# Create the sub-directories expected by run.py and symlink files known to be on AFS
mkdir samples
mkdir samples/{year}
mkdir bin
ln -s $SUBMIT_DIR/bin/eventAnalyzer bin/
ln -s $SUBMIT_DIR/libTreeAnalysis.so ./
ln -s $SUBMIT_DIR/../Producers/python/{csvname} ./

# Check if the sample is on eos
localsample=$SUBMIT_DIR/{sample_dir}/{year}/{sample}.root
abssample=$(realpath $localsample)
if ! $runninglocally && echo $abssample | grep -q ^/eos ; then
    export EOS_MGM_URL={EOS_MGM_URL}
    eos cp $abssample samples/{year}/
else
    ln -s $abssample samples/{year}/
fi

echo 'Running at:' $(date)
echo path: $(pwd)

./python/run.py -c {csvname} {runpy_args} > log.txt
runStatus=$?

echo -n $runStatus > exitStatus.txt
echo 'Done at: ' $(date) with exit status: $runStatus
gzip log.txt

echo "Files on node:"
ls -la

#delete submission scripts, so that they are not copied back (which fails sometimes)
# rm -f batchScript.sh

echo '...done at' $(date)

#note cping back is handled automatically by condor
if $runninglocally; then
  cp *.root *.txt *.gz $SUBMIT_DIR
  [ -e *.log ] && cp *log $SUBMIT_DIR
fi

exit $runStatus
'''


def get_condorsub_script( mainDir, jobFlavour='espresso' ):
    '''prepare the Condor submission script'''
    return condorsub_template.format(home=os.path.expanduser("~"), uid=os.getuid(), mainDir=mainDir, jobFlavour=jobFlavour)


def get_batch_script(sample, year, **kwargs):
    isData = sample.startswith( ('DoubleMu', 'DoubleEle', 'EGamma', 'MuEG', 'Single', 'MuonEG', 'DoubleEG', year, 'data') ) # copy-pasted from run.py
    csvname = 'samples_{year}UL_{dataMC}.csv'.format(year=year, dataMC=('Data' if isData else 'MC'))
    return script_template.format(sample=sample, year=year, csvname=csvname, **kwargs, EOS_MGM_URL=os.environ.get('EOS_MGM_URL', 'root://eosuser.cern.ch'))


def parse_args():
    parser = ArgumentParser(
        description = 'Create folders with scripts and configs to submit to HTCondor'
    )
    parser.add_argument('-A', '--analyzer'  , default='VVXAanlyzer')
    parser.add_argument('-y', '--year'      , default='2018')
    parser.add_argument('-d', '--sample-dir', default='samples')
    parser.add_argument('-j', '--flavour'   , choices=['espresso', 'microcentury', 'longlunch', 'workday', 'tomorrow', 'testmatch', 'nextweek'], default='longlunch') # TODO estimate runtime based on sample size
    parser.add_argument(      '--force'     , action='store_true', help='Delete any existing job folders')
    parser.add_argument('-o', '--output-dir', default='production', help='Base directory for the jobs')  # maybe use $(git rev-parse --short HEAD)
    parser.add_argument('--log', dest='loglevel', metavar='LEVEL', default='WARNING', help='Level for the python logging module. Can be either a mnemonic string like DEBUG, INFO or WARNING or an integer (lower means more verbose).')

    args, unknown = parser.parse_known_args()
    return args, unknown


def main():
    args, unknown_args = parse_args()
    loglevel = args.loglevel.upper() if not args.loglevel.isdigit() else int(args.loglevel)
    logging.basicConfig(format='%(levelname)s:%(module)s:%(funcName)s: %(message)s', level=loglevel)

    logging.info('args: %s', args)
    logging.info('unknown (will be passed to run.py: %s', unknown_args)

    jobDirectory = args.output_dir
    try:
        os.mkdir(jobDirectory)
    except FileExistsError as e:
        if(args.force):
            print('INFO: removing existing jobs in', jobDirectory)
            check_call(['rm', '-r', jobDirectory])
            os.mkdir(jobDirectory)
        else:
            logging.error('Job folder "%s" exists. If you want to remove and recreate its contents it run agin with --force', jobDirectory)
            raise e

    # Create condor.sub
    absJobDirectory = check_output(['realpath', jobDirectory], encoding='utf-8').strip('\n')
    with open(os.path.join(jobDirectory, 'condor.sub'), 'w') as condorsub:
        condorsub.write(get_condorsub_script( absJobDirectory, jobFlavour=args.flavour ))

    samples = os.listdir(os.path.join(args.sample_dir, args.year))
    created_jobs = 0
    for sample in ('ZZGTo4LG',):  #samples
        #Create job dir
        sample_dir = os.path.join(jobDirectory, sample)
        os.mkdir( sample_dir )

        # Create script and make it executable
        runpy_args = ' '.join(['-y', args.year, *unknown_args, '--', args.analyzer, sample])
        batchScript = get_batch_script(**vars(args), sample=sample, runpy_args=runpy_args)
        scriptPath = os.path.join(sample_dir, 'batchScript.sh')
        with open(scriptPath, 'w') as script:
            script.write(batchScript)
        os.chmod(scriptPath, os.stat(scriptPath).st_mode | stat.S_IEXEC)

        # Create logdir
        os.makedirs(os.path.join(sample_dir, 'log'), exist_ok=True)

        created_jobs += 1

    logging.info("created {:d} jobs in {:s}".format(created_jobs, absJobDirectory))


if __name__ == '__main__':
    exit(main())
