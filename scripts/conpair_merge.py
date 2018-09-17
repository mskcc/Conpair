#!/usr/bin/python

'''
@Description : This tool helps to merge concordance/contamination results from Conpair into table and generate plots in PDF format
@Created :  06/13/2018
'''

from __future__ import division
import argparse, time, os, sys, logging, re, csv, glob, subprocess

def main():
    parser = argparse.ArgumentParser(prog='conpair_merge.py', description='This tool helps to merge concordance/contamination results into table and generate plots in PDF format.', usage='%(prog)s [options]')
    parser.add_argument("-p", "--pairing", action="store", dest="inPair", required=True, type=str,
                        help="sample pairing file")
    parser.add_argument("-cl", "--cordlist", nargs='+', dest="cordList", required=True, type=str, help="A list of concordance.txt files")
    parser.add_argument('-tl', '--tamilist', nargs='+', dest="tamiList", required=True, type=str, help="A list of contamination.txt files")
    parser.add_argument('-op', '--outpre', action="store", dest="outPre", required=False, type=str, help="Prefix name for output files")
    parser.add_argument("-o", "--output", action="store", dest="outDir", required=False, type=str, help="Output directory that contains summarized tables and plots")

    args = parser.parse_args()

    # read pairing file
    normalID = []
    tumorID = []
    pair = args.inPair
    with open(pair, 'r') as pair:
        for li in pair:
            li = li.strip()
            arr_li = li.split("\t")
            normalID.append(arr_li[0])
            tumorID.append(arr_li[1])

    normalIDuniq = list(set(normalID))
    tumorIDuniq = list(set(tumorID))

    # read concordance and contamination files
    fcord = args.cordList
    ftami = args.tamiList

    # Concordance file name: s_C_LMMWH4_P001_d.Group20.rg.md.abra.printreads.pileup.s_C_04WTKC_N001_d.Group14.rg.md.abra.printreads.pileup.concordance.txt
    # Extracting information from concordance files and merge all pair-wise (N-T) concordance percentages into one table
    # output table: row - sorted normal sample ids; column - sorted tumor sample ids
    nlen, tlen = len(normalID), len(tumorID)
    cordtb = [[-1 for x in range(nlen)] for y in range(tlen)]
    for f1 in fcord:

        # assign normal and tumor ids and find their index in the 2-dimension list "cordtb"
        idx_normal = -1
        idx_tumor = -1
        for i in range(len(normalID)):
            if normalID[i] in f1:
                idx_normal = i
                break
        for j in range(len(tumorID)):
            if tumorID[j] in f1:
                idx_tumor = j
                break
        if idx_normal == -1 or -1 == idx_tumor:
            print "There is an error: Could not find sample IDs for file ",f1
            sys.exit(1)

        # read f1 file to get percentage information
        with open(f1, 'r') as cord:
            for li in cord:
                li = li.strip()
                if li.startswith("Concordance:"):
                    arr_li = li.split(" ")
                    pert = arr_li[len(arr_li) - 1]
                    pert = pert.replace("%","")
                    # assign the value of this concordance
                    cordtb[idx_normal][idx_tumor] = pert

    # get prefix name of output files
    if args.outPre:
        prefix = args.outPre
        # get output directory
        if args.outDir:
            outcord = args.outDir + "/" + prefix + "_concordance.txt"
            outcordr = args.outDir + "/" + prefix + "_concordance.R"
            outtami = args.outDir + "/" + prefix + "_contamination.txt"
            outtamir = args.outDir + "/" + prefix + "_contamination.R"
        else:
            outcord = prefix + "_concordance.txt"
            outcordr = prefix + "_concordance.R"
            outtami = prefix + "_contamination.txt"
            outtamir = prefix + "_contamination.R"
    else:
        # get output directory
        if args.outDir:
            outcord = args.outDir + "/" + "concordance.txt"
            outcordr = args.outDir + "/" + "concordance.R"
            outtami = args.outDir + "/" + "contamination.txt"
            outtamir = args.outDir + "/" + "contamination.R"
        else:
            outcord = "concordance.txt"
            outcordr = "concordance.R"
            outtami = "contamination.txt"
            outtamir = "contamination.R"
    
    # print concordance matrix (row: normal sample ID; column: tumor sample ID; each cell: percentage of concordance between two N-T samples)
    with open(outcord, 'w') as foutcord:
        # first print column names
        colname = "concordance"
        for col in tumorIDuniq:
            colname = colname + "\t" + col
        foutcord.write(colname + "\n")
        # second print each row
        counter = 0
        for rowname in normalIDuniq:
            rowstr = rowname
            row = cordtb[counter][:]
            for ele in row:
                rowstr = rowstr + "\t" + str(ele)
            foutcord.write(rowstr + "\n")
            counter += 1

    # write and run a R script to generate concordance plot in PDF format
    outcordrdir = os.path.dirname(os.path.realpath(outcordr))
    writeConcordRScript(outcordrdir, outcord, outcordr)
    cmd1 = ["Rscript", outcordr]
    execute_shell(cmd1)

    # Contamination file name: s_C_1FPU8J_N001_d.s_C_1FPU8J_P001_d.contamination.txt
    with open(outtami, 'w') as fouttami:
        # first print column names
        fouttami.write("Sample_Type" + "\t" + "Sample_ID" + "\t" + "Contamination" + "\n")

        for f2 in ftami:
            arr_f2 = f2.split("/").pop().split(".")
            tami1 = arr_f2[0]
            tami2 = arr_f2[1]
            str_normal = tami1
            str_tumor = tami2
            if tami1 in normalID:
                str_normal = tami1
            elif tami1 in tumorID:
                str_tumor = tami1
            if tami2 in normalID:
                str_normal = tami2
            elif tami2 in tumorID:
                str_tumor = tami2

            with open(f2, 'r') as tami:
                for li in tami:
                    li = li.strip()
                    arr_li = li.split(" ")
                    pert = arr_li[len(arr_li) - 1]
                    pert = pert.replace("%","")
                    if li.startswith("Normal"):
                        fouttami.write("N" + "\t" + str_normal + "\t" + pert + "\n")
                    elif li.startswith("Tumor"):
                        fouttami.write("T" + "\t" + str_tumor + "\t" + pert + "\n")

    # write and run a R script to generate contamination plot in PDF format
    outtamirdir = os.path.dirname(os.path.realpath(outtamir))
    writeContamRScript(outtamirdir, outtami, outtamir)
    cmd2 = ["Rscript", outtamir]
    execute_shell(cmd2)

def writeConcordRScript(indir, inf, outf):
    with open(outf, 'w') as rfile:
        rfile.write("library(reshape2)\n")
        rfile.write("library(ggplot2)\n")
        rfile.write("setwd(\"" + indir + "\")\n")
        rfile.write("incord <- \"" + inf + "\"\n")
        outpdf=outf.replace(".R", ".pdf")
        rfile.write("outcord <- \"" + outpdf + "\"\n")
        rfile.write("cord <- read.delim(incord, header=TRUE, row.names = 1)\n")
        rfile.write("mcord <- melt(as.matrix(cord))\n")
        rfile.write("pdf(file=outcord, width=15, height=15)\n")
        str_tmp = "ggplot(mcord, aes(Var1, Var2)) + geom_tile(aes(fill=value), colour=\"white\") + scale_fill_gradient(low=\"red\", high=\"green\") + geom_text(aes(label=value), color=\"black\", size=5) + labs(x=\"Normal Sample\", y=\"Tumor Sample\", title=\"Concordance among Smples\", fill=\"Percentage of Concordance\") + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position=\"bottom\", legend.box=\"horizontal\", plot.title=element_text(size=30, face=\"bold\", vjust=4, hjust=0.4), axis.text=element_text(size=14, face=\"bold\", colour=\"black\"), axis.title=element_text(size=20, face=\"bold\", colour=\"black\"), legend.text=element_text(size=14, face=\"bold\", colour=\"black\"), legend.title=element_text(size=20, face=\"bold\", colour=\"black\"),plot.margin=unit(c(2,2,2,2), \"lines\"))"
        rfile.write(str_tmp + "\n")
        rfile.write("dev.off()\n")

def writeContamRScript(indir, inf, outf):
    with open(outf, 'w') as rfile:
        rfile.write("library(reshape2)\n")
        rfile.write("library(ggplot2)\n")
        rfile.write("setwd(\"" + indir + "\")\n")
        rfile.write("intami <- \"" + inf + "\"\n")
        outpdf=outf.replace(".R", ".pdf")
        rfile.write("outtami <- \"" + outpdf + "\"\n")
        rfile.write("tami <- read.delim(intami, header=TRUE)\n")
        rfile.write("tami$Sample_ID <- factor(tami$Sample_ID, levels=tami$Sample_ID[order(tami$Sample_Type, -tami$Contamination)])\n")
        rfile.write("pdf(file=outtami, width=20, height=15)\n")
        str_tmp = "ggplot(tami, aes(x=Sample_ID, y=Contamination, fill=factor(Sample_Type))) + geom_bar(stat=\"identity\") + geom_hline(yintercept=2) + annotate(\"text\", min(tami$Contamination), 2, vjust=-1, hjust=-0.2, label=\"Soft Cutoff\") + geom_hline(yintercept=5) + annotate(\"text\", min(tami$Contamination), 5, vjust=-1, hjust=-0.2, label=\"Cutoff\") + labs(x=\"Samples\",y=\"Contamination(%)\",title=\"Contamination of Samples\",fill=\"Normal/Tumor\") + theme(axis.text.x=element_text(angle=45, hjust=1), legend.position=\"bottom\",legend.box=\"horizontal\", plot.title=element_text(size=30,face=\"bold\",vjust=4, hjust=0.5),axis.text=element_text(size=14,face=\"bold\",colour=\"black\"),axis.title=element_text(size=20,face=\"bold\",colour=\"black\"),legend.text=element_text(size=14,face=\"bold\",colour=\"black\"),legend.title=element_text(size=20,face=\"bold\",colour=\"black\"), plot.margin=unit(c(2,2,2,3), \"lines\")) + coord_cartesian(ylim = c(0, 100))"
        rfile.write(str_tmp + "\n")
        rfile.write("dev.off()\n")

def execute_shell(cmd):
    try:
        print >>sys.stderr, "Executing %s" % " ".join(cmd)
        subprocess.check_call(" ".join(cmd), shell=True)
    except:
        print >>sys.stderr, "Unexpected Error: %s", sys.exc_info()[0]
        sys.exit(1)

if __name__ == "__main__":
    start_time = time.time()
    main()
    end_time = time.time()
    totaltime = end_time - start_time
    print(totaltime)
    sys.exit(0)


