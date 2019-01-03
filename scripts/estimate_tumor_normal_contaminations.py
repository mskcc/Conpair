#!/usr/bin/python

# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2016) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center
# cannot be responsible for its use, misuse, or functionality.
# Version: 1.0
# Author: Ewa A Bergmann (ewa.a.bergmann@gmail.com)

# Modified on Sep-28-2018
# Purpose: Will accept list of Tumor-Normal pairs and generate PDF plot
# Author: Zuojian Tang (tangz@mskcc.org)

# Modified on Dec-10-2018
# Purpose: If Data (markers) is empty in one sample, it will not be considered for the final "table" and "pdf" files.
# Author: Zuojian Tang (tangz@mskcc.org)

# Modified on Dec-12-2018
# Purpose: If there is "pooled" sample in T/N pair, this pair will not be reported in both "table" and "pdf" files.
# Author: Zuojian Tang (tangz@mskcc.org)

# Modified on Jan-3-2019
# Purpose: If there is no considered normal/tumor samples, it will not output any summary results. This situation could happen when all normal samples are pooled normal samples.
# Purpose: If it is empty table, R codes will do nothing. In this case, there is no pdf files generated. This situtation could happen when "WARNING" message is shown in either *.concordance or *.contamination file.
# Author: Zuojian Tang (tangz@mskcc.org)

from __future__ import division
from __future__ import print_function

from collections import defaultdict
import numpy as np

import argparse, time, os, sys, subprocess
import inspect

# These modules should be in a folder under the parent of the folder containing this script
module_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))), 'modules')
# If found, add that to the PYTHONPATH, else we depend on user to define PYTHONPATH
if os.path.isdir(module_dir):
    sys.path.append(module_dir)
import MathOperations, ContaminationModel, ContaminationMarker, Genotypes

HOMOZYGOUS_P_VALUE_THRESHOLD = 0.999

def main():
    parser = argparse.ArgumentParser(prog='estimate_tumor_normal_contaminations.py', description='Program to estimate tumor-normal sample contaminations and generate summary plot in PDF format.', usage='%(prog)s [options]')
    parser.add_argument("-T", "--tumor_pileup", nargs='+', dest="t_pileups", required=True, type=str, help="A list of tumor pileup files")
    parser.add_argument("-N", "--normal_pileup", nargs='+', dest="n_pileups", required=True, type=str, help="A list of normal pileup files")
    parser.add_argument("-p", "--pairing", action="store", dest="pairing", required=True, type=str, help="sample pairing file")
    parser.add_argument("-D", "--conpair_dir", action="store", dest="conpair_dir", required=False, type=str, help="CONPAIR DIR [default: $CONPAIR_DIR]")
    parser.add_argument("-M", "--markers", action="store", dest="markers", required=False, type=str, help="MARKER FILE [default: markers for GRCh37 from $CONPAIR_DIR/data/markers/]")
    parser.add_argument("-O", "--outdir", action="store", dest="outdir", required=False, type=str, help="Output directory")
    parser.add_argument("-G", "--grid", action="store", dest="grid", required=False, type=float, default=0.01, help="GRID INTERVAL [default: 0.01]")
    parser.add_argument("-Q", "--min_mapping_quality", action="store", dest="min_mapping_quality", required=False, type=int, default=10, help="MIN MAPPING QUALITY [default: 10]")
    parser.add_argument('-op', '--outpre', action="store", dest="outPre", required=False, type=str, help="Prefix name for output files")

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

    grid_precision = args.grid
    MMQ = args.min_mapping_quality

    # get prefix name of output files
    if args.outPre:
        prefix = args.outPre
        # get output directory
        if args.outdir:
            outtami = os.path.join(args.outdir, prefix + "_contamination.txt")
            outtamir = os.path.join(args.outdir, prefix + "_contamination.R")
        else:
            outtami = prefix + "_contamination.txt"
            outtamir = prefix + "_contamination.R"
    else:
        # get output directory
        if args.outdir:
            outtami = os.path.join(args.outdir, "contamination.txt")
            outtamir = os.path.join(args.outdir, "contamination.R")
        else:
            outtami = "contamination.txt"
            outtamir = "contamination.R"

    # read pairing file
    normalID = []
    tumorID = []
    pair = args.pairing
    with open(pair, 'r') as pair:
        for li in pair:
            li = li.strip()
            arr_li = li.split("\t")
            str_normal = arr_li[0]
            str_tumor = arr_li[1]
            if "pool" in str_normal.lower() or "pool" in str_tumor.lower():
                continue
            normalID.append(str_normal)
            tumorID.append(str_tumor)

    # read pileup file lists
    ftpileup = args.t_pileups
    fnpileup = args.n_pileups

    # find pileup files for each pair in pairing file
    # and run estimate_tumor_normal_contamination
    out = []
    for i in range(len(normalID)):
        tmp_normalID = normalID[i]
        tmp_tumorID = tumorID[i]
        N_pileup = ""
        T_pileup = ""
        for fn in fnpileup:
            if tmp_normalID in fn:
                N_pileup = fn
                break
        for ft in ftpileup:
            if tmp_tumorID in ft:
                T_pileup = ft
                break
        if N_pileup == "" or T_pileup == "":
            print("There is an error: could not find pileup file for", N_pileup, "or", T_pileup)
            sys.exit(1)

        print("Estimate contamination of", T_pileup, "and", N_pileup)
        str_return = run_estimate_tumor_normal_contamination(N_pileup, T_pileup, MARKER_FILE, grid_precision, MMQ)
        arr_str_return = str_return.split("\t")
        str_ntami = arr_str_return[0]
        str_ttami = arr_str_return[1]
        arr_str_ntami = str_ntami.split(":")
        if arr_str_ntami[0] == "Normal":
            str_out = "N" + "\t" + tmp_normalID + "\t" + str(arr_str_ntami[1])
            out.append(str_out)
        arr_str_ttami = str_ttami.split(":")
        if arr_str_ttami[0] == "Tumor":
            str_out = "T" + "\t" + tmp_tumorID + "\t" + str(arr_str_ttami[1])
            out.append(str_out)

    # output into file
    '''
    Sample_Type	Sample_ID	Contamination
    N	s_C_000658_N001_d	0.692
    T	s_C_000658_T002_d	0.726
    '''
    with open(outtami, 'w') as fouttami:
        # first print column names
        fouttami.write("Sample_Type" + "\t" + "Sample_ID" + "\t" + "Contamination" + "\n")
        for tmp_item in out:
            fouttami.write(tmp_item + "\n")

    # write and run a R script to generate contamination plot in PDF format
    outtamirdir = os.path.dirname(os.path.realpath(outtamir))
    writeContamRScript(outtamirdir, outtami, outtamir)
    cmd = ["Rscript", outtamir]
    execute_shell(cmd)


def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step

def run_estimate_tumor_normal_contamination(N_pileup, T_pileup, MARKER_FILE, grid_precision, MMQ):

    Markers = ContaminationMarker.get_markers(MARKER_FILE)

    Normal_homozygous_genotype = defaultdict(lambda: defaultdict())

    checkpoints = [i for i in drange(0.0, 1.0, grid_precision)]
    checkpoints.append(0.5)
    Scores = ContaminationModel.create_conditional_likelihood_of_base_dict(checkpoints)
    checkpoints = [i for i in drange(0.0, 0.5, grid_precision)]
    checkpoints.append(0.5)

    ### PARSING THE NORMAL PILEUP FILE, CALCULATING THE LIKELIHOOD FUNCTION

    file = open(N_pileup)
    Data = []
    for line in file:
        if line.startswith("[REDUCE RESULT]"):
            continue
        pileup = ContaminationMarker.parse_mpileup_line(line, min_map_quality=MMQ)
        try:
            marker = Markers[pileup.chrom + ":" + pileup.pos]
        except:
            continue

        if pileup.Quals[marker.ref] == [] and pileup.Quals[marker.alt] == []:
            continue

        RAF = marker.RAF
        ref_basequals = pileup.Quals[marker.ref]
        alt_basequals = pileup.Quals[marker.alt]

        AA_likelihood, AB_likelihood, BB_likelihood = Genotypes.compute_genotype_likelihood(pileup.Quals[marker.ref], pileup.Quals[marker.alt], normalize=True)

        if AA_likelihood >= HOMOZYGOUS_P_VALUE_THRESHOLD:
            Normal_homozygous_genotype[pileup.chrom][pileup.pos] = {'genotype': marker.ref, 'AA_likelihood': AA_likelihood, 'AB_likelihood': AB_likelihood, 'BB_likelihood': BB_likelihood}
        elif BB_likelihood >= HOMOZYGOUS_P_VALUE_THRESHOLD:
            Normal_homozygous_genotype[pileup.chrom][pileup.pos] = {'genotype': marker.alt, 'AA_likelihood': AA_likelihood, 'AB_likelihood': AB_likelihood, 'BB_likelihood': BB_likelihood}


        p_AA, p_AB, p_BB = Genotypes.RAF2genotypeProb(RAF)
        lPAA = MathOperations.log10p(p_AA)
        lPAB = MathOperations.log10p(p_AB)
        lPBB = MathOperations.log10p(p_BB)


        priors = [lPAA*2, lPAA+lPBB, lPAA+lPAB, lPAB*2, lPAB+lPAA, lPAB+lPBB, lPBB*2, lPBB+lPAA, lPBB+lPAB]
        marker_data = [priors, ref_basequals, alt_basequals]
        Data.append(marker_data)

    file.close()
    
    # modified on Dec-10-2018
    if Data:
        D = ContaminationModel.calculate_contamination_likelihood(checkpoints, Data, Scores)
        ARGMAX = np.argmax(D)
        cont = checkpoints[ARGMAX]

        x1 = max(cont-grid_precision, 0.0)
        x2 = cont
        x3 = min(cont+grid_precision, 1.0)

        if x2 == 0.0:
            x2 += grid_precision/100
        elif x2 == 1.0:
            x2 -= grid_precision/100

        ### SEARCHING THE SPACE AROUND ARGMAX - Brent's algorithm

        optimal_val = ContaminationModel.apply_brents_algorithm(Data, Scores, x1, x2, x3)

        ### PRINTING THE NORMAL RESULTS

        str_return = "Normal:" + str(round(100.0*optimal_val, 3))
    else:
        str_return = "WARNING:" + "NA"


    ### PARSING THE TUMOR PILEUP FILE, CALCULATING THE LIKELIHOOD FUNCTION

    file = open(T_pileup)

    checkpoints = [i for i in drange(0.0, 1.0, grid_precision)]
    Data = []
    for line in file:
        if line.startswith("[REDUCE RESULT]"):
            continue
        pileup = ContaminationMarker.parse_mpileup_line(line, min_map_quality=MMQ)

        try:
            normal_hom_genotype = Normal_homozygous_genotype[pileup.chrom][pileup.pos]['genotype']
        except:
            continue

        try:
            marker = Markers[pileup.chrom + ":" + pileup.pos]
        except:
            continue

        if pileup.Quals[marker.ref] == [] and pileup.Quals[marker.alt] == []:
            continue

        RAF = marker.RAF
        ref_basequals = pileup.Quals[marker.ref]
        alt_basequals = pileup.Quals[marker.alt]

        Normal_info = Normal_homozygous_genotype[pileup.chrom][pileup.pos]
        AA_likelihood = Normal_info['AA_likelihood']
        AB_likelihood = Normal_info['AB_likelihood']
        BB_likelihood = Normal_info['BB_likelihood']
        nlPAA = MathOperations.log10p(AA_likelihood)
        nlPAB = MathOperations.log10p(AB_likelihood)
        nlPBB = MathOperations.log10p(BB_likelihood)


        p_AA, p_AB, p_BB = Genotypes.RAF2genotypeProb(RAF)
        lPAA = MathOperations.log10p(p_AA)
        lPAB = MathOperations.log10p(p_AB)
        lPBB = MathOperations.log10p(p_BB)
        priors = [lPAA+nlPAA, lPBB+nlPAA,lPAB+nlPAA, lPAB+nlPAB, lPAA+nlPAB, lPBB+nlPAB, lPBB+nlPBB, lPAA+nlPBB, lPAB+nlPBB]
        marker_data = [priors, ref_basequals, alt_basequals]
        Data.append(marker_data)

    file.close()

    # modified on Dec-10-2018
    if Data:
        D = ContaminationModel.calculate_contamination_likelihood(checkpoints, Data, Scores)
        ARGMAX = np.argmax(D)
        cont = checkpoints[ARGMAX]

        x1 = max(cont-grid_precision, 0.0)
        x2 = cont
        x3 = min(cont+grid_precision, 1.0)

        if x2 == 0.0:
            x2 += grid_precision/100
        elif x2 == 1.0:
            x2 -= grid_precision/100

        ### SEARCHING THE SPACE AROUND ARGMAX - Brent's algorithm

        optimal_val = ContaminationModel.apply_brents_algorithm(Data, Scores, x1, x2, x3)

        ### PRINTING THE TUMOR RESULTS

        str_return = str_return + "\t"+ "Tumor:" + str(round(100.0*optimal_val, 3))
    else:
        str_return = str_return + "\t"+ "WARNING:" + "NA"

    return str_return

def writeContamRScript(indir, inf, outf):

    outpdf = outf.replace(".R", ".pdf")
    with open(outf, 'w') as rfile:
        strout = """
        library(reshape2)
        library(ggplot2)
        setwd("%s")
        intami <- "%s"
        outtami <- "%s"
        tami <- read.delim(intami, header=TRUE)
        tami <- unique(tami)
        if(nrow(tami) != 0){
        tami$Sample_ID <- factor(tami$Sample_ID, levels=tami$Sample_ID[order(tami$Sample_Type, -tami$Contamination)])
        pdf(file=outtami, width=20, height=20)
        print(ggplot(tami, aes(x=Sample_ID, y=Contamination, fill=factor(Sample_Type))) + 
        geom_bar(stat="identity") +
        geom_hline(yintercept=2) + annotate("text", min(tami$Contamination), 2, vjust=-1, hjust=-0.2, label="Soft Cutoff") +
        geom_hline(yintercept=5) + annotate("text", min(tami$Contamination), 5, vjust=-1, hjust=-0.2, label="Cutoff") +
        labs(x="Samples",y="Contamination(%%)",title="Contamination of Samples", fill="Tumor/Normal ") +
        theme(plot.title=element_text(size=40,face="bold",vjust=4, hjust=0.5, margin=margin(t=10,b=20)), 
            plot.margin=unit(c(2,2,2,3), "lines"), 
            axis.title=element_text(size=28,face="bold",colour="black"), 
            axis.text=element_text(size=20,face="bold",colour="black"), 
            axis.text.x=element_text(angle=45, hjust=1, margin=margin(t=20,r=0,b=0,l=0)), 
            axis.text.y=element_text(margin=margin(t=0,r=20,b=0,l=0)), 
            legend.position="bottom", 
            legend.box="horizontal", 
            legend.spacing=unit(1,"cm"), 
            legend.text=element_text(size=24,face="bold",colour="black", margin=margin(r=50,unit="pt")),
            legend.title=element_text(size=28,face="bold",colour="black"),
            legend.key.height=unit(2,"line"), 
            legend.key.width=unit(2,"line")) +
        coord_cartesian(ylim = c(0, 100)))
        dev.off()
        }
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
