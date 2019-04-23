# Conpair

Conpair: concordance and contamination estimator for tumor–normal pairs

Conpair is a fast and robust method dedicated for human tumor-normal studies to perform concordance verification (= samples coming from the same individual), as well as cross-individual contamination level estimation in whole-genome and whole-exome sequencing experiments. Importantly, our method of estimating contamination in the tumor samples is not affected by copy number changes and is able to detect contamination levels as low as 0.1%.

* Author: Ewa A Bergmann
* Contact: ewa.a.bergmann@gmail.com
* Note: The fork at github.com/mskcc/Conpair is maintained by the Platform Informatics group at MSKCC. Please contact ckandoth@gmail.com if you have questions about it.

# Manual

**Dependencies:**

* python 2.7 or higher :     [www.python.org](https://www.python.org/)
* numpy 1.7.0 or higher :    [www.numpy.org](http://www.numpy.org/)
* scipy 0.14.0 or higher :   [www.scipy.org](http://www.scipy.org/)
* GATK 2.3 or higher :       [www.broadinstitute.org/gatk/download](http://www.broadinstitute.org/gatk/download/)
* java :                     [http://java.com](http://java.com/en/download/)
* R :                        [https://www.r-project.org](https://www.r-project.org/)

**Most common usage:**

Below is a quick guide to get started, but please use `--help` on each script to see additional options. Make sure your reference fasta has `.dict` and `.fai` index files. [Click here](http://gatkforums.broadinstitute.org/gatk/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference) to learn how.


To generate readcounts supporting ref/alt alleles at known germline sites (GATK required):
```
python scripts/run_gatk_pileup_for_sample.py -B TUMOR.bam -O TUMOR.pileup -R /path/to/reference.fa -G /path/to/gatk.jar -M data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt
```

To verify concordance between one or more pairs of tumor/normal samples, and generate a PDF summary:
```
python scripts/verify_concordances.py -T TUMOR.pileup -N NORMAL.pileup -p pairing.txt -M data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt -O output_dir
```

To estimate contamination levels in each tumor and normal:

```
python scripts/estimate_tumor_normal_contaminations.py -T TUMOR.pileup -N NORMAL.pileup -p pairing.txt -M data/markers/GRCh37.autosomes.phase3_shapeit2_mvncall_integrated.20130502.SNV.genotype.sselect_v4_MAF_0.4_LD_0.8.txt -O output_dir
```

# Output files

Sample outputs can be found under `data/example` in this repo.

# Interpretation

**Concordance**

To eliminate the effect of copy number variation on the concordance levels, we recommend using the -H (--normal_homozygous_markers_only) flag.
If two samples are concordant the expected concordance level should be close to 99-100%.
For discordant samples concordance level should be close to 40%.
You can observe slighly lower concordance (80-99%) in presence of contamination and/or copy number changes (if the -H option wasn't used) in at least one of the samples.

**Contamination**

Even a low contamination level (such as 0.5%) in the tumor sample will effect somatic mutation calls, resulting in additional false-positives at variant allele fractions around that contamination level. Cross-individual contamination in the normal sample is less likely to affect somatic mutation callers, but may affect allele-specific somatic copy-number callers.
