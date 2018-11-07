#!/usr/bin/python

'''
@Description : This tool helps to merge concordance/contamination results from Conpair into table and generate plots in PDF format
@Created :  06/13/2018
@modified : 09/17/2018
purpose : There is a possibility that one normal sample matches to multiple tumor samples. The code has been modified to fit this situation.
@modified : 11/05/2018
@purpose : If normal and tumor samples do not have shared markers, there will be "WARNING" message in *.concordance.txt file. This pair will not be reported in summarized pdf file.
@purpose : If it is either pooled normal or tumor sample, this pair will be reported in summarized pdf file.
@modified : 11/07/2018
@purpose : If either normal or tumor has lower coverage, it will have "WARNING" message in *.contamination.txt file. This pair will be reported in summarized pdf file
'''

from __future__ import division
import argparse, time, os, sys, logging, re, csv, glob, subprocess
import inspect

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
    normalIDuniq = []
    tumorIDuniq = []
    pair = args.inPair
    with open(pair, 'r') as pair:
        for li in pair:
            li = li.strip()
            arr_li = li.split("\t")
            tmp_normalID = arr_li[0]
            tmp_tumorID = arr_li[1]
            if "pool" not in tmp_normalID.lower() and tmp_normalID not in normalID:
                normalIDuniq.append(tmp_normalID)
            if "pool" not in tmp_tumorID.lower() and tmp_tumorID not in tumorID:
                tumorIDuniq.append(tmp_tumorID)
            normalID.append(tmp_normalID)
            tumorID.append(tmp_tumorID)

    # read concordance and contamination files
    fcord = args.cordList
    ftami = args.tamiList

    # Concordance file name: s_C_LMMWH4_P001_d.Group20.rg.md.abra.printreads.pileup.s_C_04WTKC_N001_d.Group14.rg.md.abra.printreads.pileup.concordance.txt
    # Extracting information from concordance files and merge all pair-wise (N-T) concordance percentages into one table
    # output table: row - sorted normal sample ids; column - sorted tumor sample ids
    nlen, tlen = len(normalIDuniq), len(tumorIDuniq)
    cordtb = [[-1 for y in range(tlen)] for x in range(nlen)]
    for f1 in fcord:
        # assign normal and tumor ids and find their index in the 2-dimension list "cordtb"
        idx_normal = -1
        idx_tumor = -1
        for i in range(len(normalIDuniq)):
            if normalIDuniq[i] in f1:
                idx_normal = i
                break
        for j in range(len(tumorIDuniq)):
            if tumorIDuniq[j] in f1:
                idx_tumor = j
                break
        if idx_normal == -1 or idx_tumor == -1:
            if "pool" not in f1.lower():
                print "There is an error: Could not find sample IDs in the pairing file ", f1
                sys.exit(1)
            else:
                continue
        # read f1 file to get percentage information
        with open(f1, 'r') as cord:
            for li in cord:
                li = li.strip()
                if li.startswith("WARNING:"):
                    cordtb[idx_normal][idx_tumor] = "NA"
                    break
                elif li.startswith("Concordance:"):
                    arr_li = li.split(" ")
                    pert = arr_li[len(arr_li) - 1]
                    pert = pert.replace("%","")
                    # assign the value of this concordance
                    cordtb[idx_normal][idx_tumor] = pert
                    break

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
            colname = colname + "\t" + str(col)
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

            if "pool" in str_normal.lower() or "pool" in str_tumor.lower():
                continue
            str_out = ""
            with open(f2, 'r') as tami:
                for li in tami:
                    li = li.strip()
                    arr_li = li.split(" ")
                    pert = arr_li[len(arr_li) - 1]
                    pert = pert.replace("%","")
                    if li.startswith("WARNING:"):
                        break
                    if li.startswith("Normal"):
                        str_out = "N" + "\t" + str_normal + "\t" + pert + "\n"
                        # fouttami.write("N" + "\t" + str_normal + "\t" + pert + "\n")
                    elif li.startswith("Tumor"):
                        str_out = str_out + "T" + "\t" + str_tumor + "\t" + pert + "\n"
                        # fouttami.write("T" + "\t" + str_tumor + "\t" + pert + "\n")
                fouttami.write(str_out)

    # write and run a R script to generate contamination plot in PDF format
    outtamirdir = os.path.dirname(os.path.realpath(outtamir))
    writeContamRScript(outtamirdir, outtami, outtamir)
    cmd2 = ["Rscript", outtamir]
    execute_shell(cmd2)

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
        tami$Sample_ID <- factor(tami$Sample_ID, levels=tami$Sample_ID[order(tami$Sample_Type, -tami$Contamination)])
        pdf(file=outtami, width=20, height=20)
        ggplot(tami, aes(x=Sample_ID, y=Contamination, fill=factor(Sample_Type))) + 
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
        coord_cartesian(ylim = c(0, 100))
        dev.off()
        """
        strout = inspect.cleandoc(strout)
        strout = strout % (indir, inf, outpdf)
        rfile.write(strout)

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


