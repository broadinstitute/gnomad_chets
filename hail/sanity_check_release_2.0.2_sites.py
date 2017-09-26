from hail import *
from utils import *
from sites_vcf import run_sites_sanity_checks
import time
import argparse


#Run sanity checks on split vdses?

def main(args):
    hc = HailContext()

    pops = GENOME_POPS if args.genomes else EXOME_POPS
    file_path = 'gs://gnomad/release/2.0.2/vds/genomes/gnomad.genomes.r2.0.2.sites.vds' if args.genomes else 'gs://gnomad/release/2.0.2/vds/exomes/gnomad.exomes.r2.0.2.sites.vds'

    vds = hc.read(file_path)
    vds = vds.repartition(2000, shuffle=False)
    sanity_check = run_sites_sanity_checks(vds, pops)
    if args.slack_channel: send_snippet(args.slack_channel, sanity_check, 'sanity_%s.txt' % time.strftime("%Y-%m-%d_%H:%M"))



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--exomes', help='Input VDS is exomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--genomes', help='Input VDS is genomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    args = parser.parse_args()

    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified')

    main(args)