#!/usr/bin/env python2.7

# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2016) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 0.15
# Author: Ewa A Bergmann (ewa.a.bergmann@gmail.com)

# Modified on Oct-08-2018
# Purpose: Will accept list of Tumor-Normal pairs and generate PDF plot
# Author: Zuojian Tang (tangz@mskcc.org)

from __future__ import division
from __future__ import print_function

import sys
import os
import optparse
import math
from collections import defaultdict

import argparse, time, subprocess
import inspect

# These modules should be in a folder under the parent of the folder containing this script
module_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'modules')
# If found, add that to the PYTHONPATH, else we depend on user to define PYTHONPATH
if os.path.isdir(module_dir):
    sys.path.append(module_dir)
from ContaminationModel import *
from ContaminationMarker import *

def main():
    parser = argparse.ArgumentParser(prog='verify_concordances.py',
                                     description='Program to verify tumor-normal sample concordances and generate summary plot in PDF format.',
                                     usage='%(prog)s [options]')
    parser.add_argument("-T", "--tumor_pileup", nargs='+', dest="t_pileups", required=True, type=str,
                        help="A list of tumor pileup files")
    parser.add_argument("-N", "--normal_pileup", nargs='+', dest="n_pileups", required=True, type=str,
                        help="A list of normal pileup files")
    parser.add_argument("-p", "--pairing", action="store", dest="pairing", required=True, type=str,
                        help="sample pairing file")
    parser.add_argument("-D", "--conpair_dir", action="store", dest="conpair_dir", required=False, type=str,
                        help="CONPAIR DIR [default: $CONPAIR_DIR]")
    parser.add_argument("-M", "--markers", action="store", dest="markers", required=False, type=str,
                        help="MARKER FILE [default: markers for GRCh37 from $CONPAIR_DIR/data/markers/]")
    parser.add_argument("-O", "--outdir", action="store", dest="outdir", required=False, type=str,
                        help="Output directory")
    parser.add_argument("-C", "--min_cov", action="store", dest="min_cov", required=False, type=int, default=10,
                        help="MIN COVERAGE TO CALL GENOTYPE [default: 10]")
    parser.add_argument("-Q", "--min_mapping_quality", action="store", dest="min_mapping_quality", required=False,
                        type=int, default=10, help="MIN MAPPING QUALITY [default: 10]")
    parser.add_argument("-B", "--min_base_quality", action="store", dest="min_base_quality", required=False,
                        type=int, default=20, help="MIN BASE QUALITY [default: 20]")
    parser.add_argument("-H", "--normal_homozygous_markers_only", action="store_true", required=False,
                        help="USE ONLY MARKERS THAT ARE HOMOZYGOUS IN THE NORMAL SAMPLE (concordance will not be affected by CNV)")
    parser.add_argument('-op', '--outpre', action="store", dest="outPre", required=False, type=str,
                        help="Prefix name for output files")

    args = parser.parse_args()

    if args.conpair_dir:
        CONPAIR_DIR = args.conpair_dir
    elif 'CONPAIR_DIR' in os.environ:
        CONPAIR_DIR = os.environ['CONPAIR_DIR']
    else:
        CONPAIR_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

    if args.markers:
        MARKER_FILE = args.markers
    else:
        MARKER_FILE = os.path.join(CONPAIR_DIR, 'data', 'markers', 'GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt')

    COVERAGE_THRESHOLD = args.min_cov
    MMQ = args.min_mapping_quality
    MBQ = args.min_base_quality
    AA_BB_only = args.normal_homozygous_markers_only

    # get prefix name of output files
    if args.outPre:
        prefix = args.outPre
        # get output directory
        if args.outdir:
            outcord = os.path.join(args.outdir, prefix + "_concordance.txt")
            outcordr = os.path.join(args.outdir, prefix + "_concordance.R")
        else:
            outcord = prefix + "_concordance.txt"
            outcordr = prefix + "_concordance.R"
    else:
        # get output directory
        if args.outdir:
            outcord = os.path.join(args.outdir, "concordance.txt")
            outcordr = os.path.join(args.outdir, "concordance.R")
        else:
            outcord = "concordance.txt"
            outcordr = "concordance.R"

    # read pairing file
    normalID = []
    tumorID = []
    normalIDuniq = []
    tumorIDuniq = []
    pair = args.pairing
    with open(pair, 'r') as pair:
        for li in pair:
            li = li.strip()
            arr_li = li.split("\t")
            tmp_normalID = arr_li[0]
            tmp_tumorID = arr_li[1]
            if tmp_normalID not in normalID:
                normalIDuniq.append(tmp_normalID)
            if tmp_tumorID not in tumorID:
                tumorIDuniq.append(tmp_tumorID)
            normalID.append(tmp_normalID)
            tumorID.append(tmp_tumorID)

    # read pileup file lists
    ftpileup = args.t_pileups
    fnpileup = args.n_pileups

    # run run_verify_concordance for each pair of T/N pair
    nlen, tlen = len(normalIDuniq), len(tumorIDuniq)
    total = nlen*tlen
    print(nlen, "vs", tlen, "and total=", total)
    count = 0
    cordtb = [[-1 for t in range(tlen)] for n in range(nlen)]
    for row in range(len(normalIDuniq)):
        nid = normalIDuniq[row]
        N_pileup = ""
        # search for normal pileup file
        for npileup in fnpileup:
            npath, nname = os.path.split(npileup)
            if nname.startswith(nid):
                N_pileup = npileup
                break
        # for each tumor pileup file, do concordance
        for col in range(len(tumorIDuniq)):
            tid = tumorIDuniq[col]
            T_pileup = ""
            # search for tumor pileup file
            for tpileup in ftpileup:
                tpath, tname = os.path.split(tpileup)
                if tname.startswith(tid):
                    T_pileup = tpileup
                    break
            if N_pileup == "" or T_pileup == "":
                print ("There is an error: Could not find pileup file. (", nid, tid, ")")
                sys.exit(1)
            count += 1
            print(count,"/",total)
            str_return = run_verify_concordance(N_pileup, T_pileup, MARKER_FILE, MMQ, MBQ, COVERAGE_THRESHOLD, AA_BB_only)
            arr_str_return = str_return.split("\t")
            str_cord = arr_str_return[0]
            arr_cord = str_cord.split(":")
            if arr_cord[0] == "Concordance":
                pert = arr_cord[1]
                pert = pert.replace("%", "")
                # assign the value of this concordance
                cordtb[row][col] = pert

    # print concordance matrix (row: normal sample ID; column: tumor sample ID; each cell: percentage of concordance between two N-T samples)
    with open(outcord, 'w') as foutcord:
        # first print column names
        colname = "concordance"
        for col in tumorIDuniq:
            colname = colname + "\t" + str(col)
        foutcord.write(colname + "\n")
        # second print each row
        for row_num in range(len(normalIDuniq)):
            rowname = normalIDuniq[row_num]
            rowstr = rowname
            row = cordtb[row_num][:]
            for ele in row:
                rowstr = rowstr + "\t" + str(ele)
            foutcord.write(rowstr + "\n")

    # write and run a R script to generate concordance plot in PDF format
    outcordrdir = os.path.dirname(os.path.realpath(outcordr))
    writeConcordRScript(outcordrdir, outcord, outcordr)
    cmd = ["Rscript", outcordr]
    execute_shell(cmd)

def run_verify_concordance (N_pileup, T_pileup, MARKER_FILE, MMQ, MBQ, COVERAGE_THRESHOLD, AA_BB_only):

    Markers = get_markers(MARKER_FILE)

    Normal_genotype_likelihoods = genotype_likelihoods_for_markers(Markers, N_pileup, min_map_quality=MMQ, min_base_quality=MBQ)
    Tumor_genotype_likelihoods = genotype_likelihoods_for_markers(Markers, T_pileup, min_map_quality=MMQ, min_base_quality=MBQ)

    concordant = 0
    discordant = 0
    for m in Markers:
        NL = Normal_genotype_likelihoods[m]
        TL = Tumor_genotype_likelihoods[m]
        if NL is None or TL is None:
            continue
        if NL['coverage'] < COVERAGE_THRESHOLD or TL['coverage'] < COVERAGE_THRESHOLD:
            continue
        if AA_BB_only:
            if NL['likelihoods'].index(max(NL['likelihoods'])) == 1:
                continue
        if NL['likelihoods'].index(max(NL['likelihoods'])) == TL['likelihoods'].index(max(TL['likelihoods'])):
            concordant += 1
        else:
            discordant += 1

    if concordant+discordant == 0:
        print('WARNING: There are no shared markers between the tumor and the normal samples that meet the specified coverage requirements ({0})\nIs the coverage of your samples high enough?\nExiting...'.format(COVERAGE_THRESHOLD))
        sys.exit(0)

    str_return = "Concordance:"+ str(round(100.0*float(concordant)/(concordant+discordant), 2)) + "%"
    str_return = str_return + "\t" + \
                 "Based on " + str(concordant+discordant) + "/" + str(len(Markers)) + " markers (coverage per marker threshold:" + str(COVERAGE_THRESHOLD) + " reads)" + "\t" + \
                 "Minimum mappinq quality:" + str(MMQ) + "\t" + \
                 "Minimum base quality:" + str(MBQ)
    return str_return

def writeConcordRScript(indir, inf, outf):

    outpdf = outf.replace(".R", ".pdf")
    with open(outf, 'w') as rfile:
        strout = """
        library(reshape2)
        library(ggplot2)
        setwd("%s")
        incord <- "%s"
        outcord <- "%s"
        cord <- read.delim(incord, header=TRUE, row.names = 1)
        mcord <- melt(as.matrix(cord))
        pdf(file=outcord, width=20, height=20)
        ggplot(mcord, aes(Var1, Var2)) + 
        geom_tile(aes(fill=value), colour="white") + 
        scale_fill_gradient(low="red", high="green", limits=c(0, 100)) + 
        geom_text(aes(label=value), color="black", size=7) +
        labs(x="\\nNormal Sample\\n", y="Tumor Sample\\n", title="Concordance among Samples") + 
        theme(plot.title=element_text(size=40, face="bold", vjust=4, hjust=0.4, margin=margin(t=10,b=20)), 
            plot.margin=unit(c(2,2,2,2), "lines"), 
            axis.title=element_text(size=28, face="bold", colour="black"), 
            axis.text=element_text(size=20, face="bold", colour="black"), 
            axis.text.x=element_text(angle=45, hjust=1, margin=margin(t=20,r=0,b=0,l=0)), 
            axis.text.y=element_text(margin=margin(t=0,r=20,b=0,l=0)), 
            legend.position="bottom", 
            legend.box="horizontal", 
            legend.spacing=unit(1,"cm"), 
            legend.text=element_text(size=24, face="bold", colour="black", margin=margin(t=10)),
            legend.title=element_text(size=28, face="bold", colour="black")) + 
        guides(fill=guide_colorbar(title="Percentage of Concordance",
            label.position="bottom", 
            title.position="left", 
            title.vjust=1, 
            title.hjust=500, 
            label=TRUE, 
            draw.ulim=TRUE, 
            draw.llim=TRUE, 
            frame.colour="black", 
            ticks=TRUE, 
            barwidth=15, 
            barheight=1.5, 
            direction='horizontal'))
        dev.off()
        """
        strout = inspect.cleandoc(strout)
        strout = strout % (indir, inf, outpdf)
        rfile.write(strout)

def execute_shell(cmd):
    try:
        # print >>sys.stderr, "Executing %s" % " ".join(cmd)
        print("Executing %s" % " ".join(cmd), file=sys.stderr)
        subprocess.check_call(" ".join(cmd), shell=True)
    except:
        # print >>sys.stderr, "Unexpected Error: %s", sys.exc_info()[0]
        print("Unexpected Error: %s" % sys.exc_info()[0], file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    print(totaltime)
    sys.exit(0)
                  

