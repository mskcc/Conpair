export SHELL:=/bin/bash
.ONESHELL:
UNAME:=$(shell uname)

define help
endef
export help
help:
	@printf "$$help"
.PHONY : help

# ~~~~~ Install Dependencies ~~~~~ #
# need to use Python 2.7 because Python 3 gives different results
export PATH:=$(CURDIR)/conda/bin:$(CURDIR)/bin:$(PATH)
unexport PYTHONPATH
unexport PYTHONHOME

ifeq ($(UNAME), Darwin)
CONDASH:=Miniconda2-4.7.12.1-MacOSX-x86_64.sh
endif

ifeq ($(UNAME), Linux)
CONDASH:=Miniconda2-4.7.12.1-Linux-x86_64.sh
endif

CONDAURL:=https://repo.anaconda.com/miniconda/$(CONDASH)

conda:
	echo ">>> Setting up conda..."
	wget "$(CONDAURL)"
	bash "$(CONDASH)" -b -p conda
	rm -f "$(CONDASH)"

# https://github.com/mskcc/roslin-variant/blob/2.6.x/build/containers/conpair/0.3.3/Dockerfile
install: conda
	pip install scipy==1.1.0 numpy==1.15.4

# python ../Conpair/scripts/verify_concordances.py -p pairing.txt -N ... -T ...
# python ../Conpair/scripts/estimate_tumor_normal_contaminations.py -p pairing.txt -N ... -T ...
# export GATK_JAR=/path/to/gatk.jar
TUMOR_PILEUP:=$(shell head -1 tumor_pileups.1.txt)
NORMAL_PILEUP:=$(shell head -1 normal_pileups.1.txt)
test-run1:
	export CONPAIR_DIR=$(CURDIR)
	python2.7 scripts/verify_concordance.py \
	--tumor_pileup $(TUMOR_PILEUP) \
	--normal_pileup $(NORMAL_PILEUP) && \
	python2.7 scripts/estimate_tumor_normal_contamination.py \
	--tumor_pileup $(TUMOR_PILEUP) \
	--normal_pileup $(NORMAL_PILEUP)
# 0.782
# Based on 188/7387 markers (coverage per marker threshold: 10 reads)
# Minimum mappinq quality: 10
# Minimum base quality: 20
# Normal sample contamination level: 0.113%
# Tumor sample contamination level: 0.0%
THREADS:=8
NUM_TUMORS:=1
NUM_NORMALS:=5
run:
	export CONPAIR_DIR=$(CURDIR)
	python run.py "$(THREADS)" "$(NUM_TUMORS)" "$(NUM_NORMALS)"

bash:
	bash