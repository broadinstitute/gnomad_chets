#!/usr/bin/env python
import argparse
import sys
import subprocess
import os
import tempfile

standard_scripts = ['~/pending_research/exac/gnomad_qc/hail/utils.py',
                    '~/pending_research/exac/gnomad_qc/hail/variantqc.py',
                    '~/pending_research/exac/gnomad_qc/hail/resources.py',
                    '~/pending_research/exac/gnomad_qc/hail/slack_creds.py',
                    '~/pending_research/exac/gnomad_qc/hail/slack_utils.py']
standard_scripts = [os.path.expanduser(x) for x in standard_scripts]


def main(args, pass_through_args):
    temp_py = None
    if args.script is None:
        if args.inline is None:
            print >> sys.stderr, 'Either --script or --inline is required. Exiting.'
            sys.exit(1)
        if 'print' not in args.inline:
            pass
        temp_py = tempfile.mkstemp(suffix='.py')
        with open(temp_py, 'w') as temp_py_f:
            temp_py_f.write(args.inline)
        script = temp_py
    else:
        script = args.script

    print >> sys.stderr, 'Running %s on %s' % (args.script, args.cluster)

    hash_string = ''
    try:
        with open(os.path.expanduser('~/Dropbox/src/misc/hail_latest_hash.txt')) as f:
            hash_string = f.read().strip()
    except IOError:
        pass

    if not hash_string:
        subprocess.check_output(['gsutil', 'cat', 'gs://hail-common/latest-hash.txt'])

    if not hash_string:
        print >> sys.stderr, 'Could not get hash string'
        sys.exit(1)

    if args.jar is not None:
        jar = args.jar
        jar_file = os.path.basename(jar)
    else:
        jar_file = 'hail-hail-is-master-all-spark2.0.2-%s.jar' % hash_string
        jar = 'gs://hail-common/%s' % jar_file

    pyfiles = [args.pyhail if args.pyhail is not None else 'gs://hail-common/pyhail-hail-is-master-%s.zip' % hash_string]

    if args.add_scripts:
        pyfiles.extend([os.path.expanduser(x) for x in args.add_scripts.split(',')])
    pyfiles.extend(standard_scripts)

    print >> sys.stderr, 'Using JAR: %s and files:\n%s' % (jar, '\n'.join(pyfiles))

    job = ['gcloud', 'dataproc', 'jobs', 'submit', 'pyspark', script,
           '--cluster', args.cluster,
           '--files=%s' % jar,
           '--py-files=%s' % ','.join(pyfiles),
           '--properties=spark.driver.extraClassPath=./%s,spark.executor.extraClassPath=./%s' % (jar_file, jar_file),
           '--driver-log-levels', 'root=FATAL,is.hail=INFO'
    ]
    if pass_through_args is not None:
        job.append('--')
        job.extend(pass_through_args)

    subprocess.check_output(job)
    if temp_py is not None:
        os.remove(temp_py)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--script', '--input', '-i', help='Script to run')
    parser.add_argument('--inline', help='Inline script to run')
    parser.add_argument('--cluster', help='Which cluster to run on', default='exomes')
    parser.add_argument('--jar', help='Jar file to use')
    parser.add_argument('--pyhail', help='Pyhail zip file to use')
    parser.add_argument('--add_scripts', help='Comma-separated list of additional scripts to add')
    parser.add_argument('--skip_other_scripts', help='Do not upload helper scripts (utils.py, etc)', action='store_true')
    args, pass_through_args = parser.parse_known_args()
    main(args, pass_through_args)
