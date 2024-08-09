#!/bin/env python3

# Creates HTCondor jobs for an EventAnayzer

import os
import stat
import re
import logging
from subprocess import check_call, check_output #, CalledProcessError
from argparse import ArgumentParser


thresholds = [(1.2e9, 'espresso'), (5e9, 'microcentury'), (10e9, 'longlunch'), (50e9, 'workday')]  # 50GB is the max file size on EOS

condorsub_template = '''\
executable              = $(directory)/batchScript.sh
arguments               = {mainDir}/$(directory) $(ClusterId)$(ProcId)
output                  = log/$(ClusterId).$(ProcId).out
error                   = log/$(ClusterId).$(ProcId).err
log                     = log/$(ClusterId).$(ProcId).log
Initialdir              = $(directory)
environment             = "CMSSW_BASE={CMSSW_BASE}"
{requirements}
request_memory          = 4000M
+JobFlavour             = "{jobFlavour}"

x509userproxy           = {home}/x509up_u{uid}

#https://www-auth.cs.wisc.edu/lists/htcondor-users/2010-September/msg00009.shtml
periodic_remove         = JobStatus == 5

ShouldTransferFiles     = YES
'''

script_template = '''\
#!/bin/bash
set -ex
set -o pipefail

if [ -z ${{_CONDOR_SCRATCH_DIR+x}} ]; then
  runninglocally=true
  _CONDOR_SCRATCH_DIR=$(mktemp -d)
  SUBMIT_DIR="$(pwd)"
else
  runninglocally=false
  SUBMIT_DIR="$1"
fi

treeanalysis_dir={CMSSW_BASE}/src/VVXAnalysis/TreeAnalysis

printf "Access to AFS: " ; [ -e /afs ]      && echo Y || echo N
printf "Access to EOS: " ; [ -e /eos/user ] && echo Y || echo N
printf "Access to treeanalyis_dir (%s): " "$treeanalysis_dir" ; [ -e "$treeanalysis_dir" ] && echo Y || echo N

source /cvmfs/cms.cern.ch/cmsset_default.sh || echo "WARN: cmsset_default.sh failed with status $?"
cd "$treeanalysis_dir"
cmsenv
cd "$_CONDOR_SCRATCH_DIR"

# Create the sub-directories expected by run.py and symlink files known to be on AFS
ln -s $treeanalysis_dir/data

echo 'Running at:' $(date)
echo path: $(pwd)
echo "OS            :" $(uname -a)

runStatus=0
$treeanalysis_dir/python/run.py -e -d "$treeanalysis_dir"/{samplesdir} {runpy_extra_args} -- {runpy_args} &>run.log || runStatus=$?

echo -n $runStatus > exitStatus.txt
echo 'run.py done at: %s, with exit status: %d' "$(date)" $runStatus

echo "Files on node:"
ls -la

#delete submission scripts, so that they are not copied back (which fails sometimes)
rm batchScript.sh || true

#note cping back is handled automatically by condor
if $runninglocally; then
  cp -r results *.out *.err *.log exitStatus.txt $SUBMIT_DIR
  cd $SUBMIT_DIR
  # rm -r $_CONDOR_SCRATCH_DIR
else
  cp -r results $SUBMIT_DIR/ || echo "ERROR: copying back results failed"
  cp *.out *.err *.log $SUBMIT_DIR/  || echo "ERROR: copying back the out and err files failed"
fi

echo 'Batch script done at' $(date)

exit $runStatus
'''


def get_os_requirements():
    # Use the same Red Hat release as the submission machine
    # This code is specific to rhel and will need updating if the OS changes (or if another version comes out)
    with open('/etc/redhat-release') as f:
        release = f.read().rstrip('\n')
    logging.debug('OS release on this machine: %s', release)
    if "release 7" in release:
        req = 'MY.SingularityImage     = "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-cat/cmssw-lxplus/cmssw-el7-lxplus:latest/"'
        # req = "requirements = (OpSysAndVer =?= \"CentOS7\")"
    elif "release 8" in release: #use a Singularity container
        req = "MY.WantOS               = \"el8\""
    elif "release 9" in release:
        req = "requirements            = (OpSysAndVer =?= \"AlmaLinux9\")"
    else:
        raise RuntimeError('Unknown Red Hat release "%s"' %(release))
    logging.debug('therefore this job has: %s', req)
    return req


def get_condorsub_script(mainDir, jobFlavour='espresso', requirements=''):
    '''prepare the Condor submission script'''

    return condorsub_template.format(home=os.path.expanduser("~")
                                     , uid=os.getuid()
                                     , mainDir=mainDir
                                     , jobFlavour=jobFlavour
                                     , requirements=requirements
                                     , CMSSW_BASE=os.environ['CMSSW_BASE']
    )


def get_batch_script(sample, year, samples_dir=None, **kwargs):
    isData = sample.startswith( ('DoubleMu', 'DoubleEle', 'EGamma', 'MuEG', 'Single', 'MuonEG', 'DoubleEG', year, 'data') ) # copy-pasted from run.py
    dataMC = 'Data' if isData else 'MC'
    if(samples_dir is None):
        samples_dir = "samples/"+dataMC
    return script_template.format(
        sample=sample
        , year=year
        , samplesdir=samples_dir
        , **kwargs
        , EOS_MGM_URL=os.environ.get('EOS_MGM_URL', 'root://eosuser.cern.ch')
        , CMSSW_BASE=os.environ['CMSSW_BASE']
    )


def parse_args():
    parser = ArgumentParser(
        description = 'Create job folders to submit the EventAnalyzer step to HTCondor'
        , epilog='Unknown arguments will be forwarded to run.py in the batch scripts'
    )
    parser.add_argument('-A', '--analyzer'  , default='VVXAnalyzer', help='Default: %(default)s')
    parser.add_argument('-y', '--year'      , default='2018'       , help='Default: %(default)s')
    parser.add_argument('-d', '--samples-dir', default='samples'    , help='Sample location, similar to run.py (default: %(default)s)')
    parser.add_argument('-r', '--regions'   , default='SR4P,CR3P1F,CR2P2F,SR3P,CR110,CR101,CR011,CR100,CR010,CR001,CR000,SR2P,CRLFR', help='Passed verbatim to run.py. Default: %(default)s')
    parser.add_argument(      '--force'     , action='store_true'  , help='Delete any existing job folders')
    parser.add_argument('-o', '--output-dir', default='production' , help='Base directory for the jobs (default: %(default)s)')  # maybe use $(git rev-parse --short HEAD)
    parser.add_argument('--log', dest='loglevel', metavar='LEVEL', default='WARNING', help='Level for the python logging module. Can be either a mnemonic string like DEBUG, INFO or WARNING or an integer (lower means more verbose).')

    args, unknown = parser.parse_known_args()
    return args, unknown


def main(args, unknown_args):
    logging.info('args: %s', args)
    logging.info('unknown options (will be passed to run.py): %s', unknown_args)

    os_requirements = get_os_requirements()

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    samples = [s.rstrip('.root') for s in os.listdir(os.path.join(args.samples_dir, args.year)) if s.endswith('.root')]
    created_jobs = 0
    for sample in samples: #('ZZGTo4LG',): #
        regex_part_match = re.search('_part(\d+)(of\d)?$', sample)
        if regex_part_match:
            part_number = int(regex_part_match.group(1))
            # "Chunk" is used so that production management scripts (e.g. haddChunks.py) will also work here
            sample_base = sample.split('_part')[0]
            job_name    = '%s_Chunk%d' %(sample_base, part_number)
        else:
            sample_base = sample
            job_name    = sample

        sample_dir = os.path.join(output_dir, args.year, sample_base)
        job_dir    = os.path.join(sample_dir, job_name)

        # Create sample dir, if necessary
        os.makedirs(sample_dir, exist_ok=True)

        # Create the dir for this specific job (=chunk)
        if(args.force):
            if(os.path.exists(job_dir)):
                check_call(['rm', '-r', job_dir])
        os.mkdir(job_dir)

        # Create condor.sub
        condorsub_path = os.path.join(sample_dir, 'condor.sub')
        if(not os.path.exists(condorsub_path)):
            # For multipart files there is a common condor sub, created for the first part
            abs_job_dir = check_output(['realpath', sample_dir], encoding='utf-8').strip('\n')
            flavour = get_job_flavour(os.path.join(args.samples_dir, args.year, sample+'.root'))
            condorsub_script = get_condorsub_script(abs_job_dir, jobFlavour=flavour, requirements=os_requirements)
            with open(condorsub_path, 'w') as condorsub:
                condorsub.write(condorsub_script)

        # Create script and make it executable
        runpy_args = ' '.join([args.analyzer, sample])
        runpy_extra_args = ' '.join(['-y', args.year, '-r', args.regions, *unknown_args])
        batchScript = get_batch_script(**vars(args), sample=sample, runpy_args=runpy_args, runpy_extra_args=runpy_extra_args)
        scriptPath = os.path.join(job_dir, 'batchScript.sh')
        with open(scriptPath, 'w') as script:
            script.write(batchScript)
        os.chmod(scriptPath, os.stat(scriptPath).st_mode | stat.S_IEXEC)

        # Create logdir
        os.mkdir(os.path.join(job_dir, 'log'))#, exist_ok=True)

        created_jobs += 1

    logging.info("created {:d} jobs in {:s}".format(created_jobs, output_dir))


def get_job_flavour(sample):
    size = os.stat(sample).st_size
    for thr, flav in thresholds:
        if(size < thr):
            flavour = flav
            break
    else:
        raise RuntimeError('sample %s has size greater than max (%d)' %(sample, thr))

    return flavour


if __name__ == '__main__':
    args, unknown_args = parse_args()
    loglevel = args.loglevel.upper() if not args.loglevel.isdigit() else int(args.loglevel)
    logging.basicConfig(format='%(levelname)s:%(module)s:%(funcName)s: %(message)s', level=loglevel)

    exit(main(args, unknown_args))
