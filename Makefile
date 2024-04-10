ROOT=/net/topmed11/working/porchard/rnaseq-make-scan-input
DATA=$(ROOT)/data
WORK=$(ROOT)/work
BIN=$(ROOT)/bin

ANALYSIS=$(WORK)/$@

.PHONY: all

define NL


endef

##### DATA #####
data: gtf unrelated genotype-pca metadata pass-qc gene-expression mappability

gtf:
	mkdir -p $(DATA)/$@ && cd $(DATA)/$@ && wget https://personal.broadinstitute.org/francois/topmed/gencode.v30.GRCh38.ERCC.genes.collapsed_only.gtf.gz

# ancestry:
# 	mkdir -p $(DATA)/$@
# 	cp /net/topmed8/working/jomiwang/freeze9b/LAI/releasing/topmed.global.ancestry.txt $(DATA)/$@/
# 	printf "NWD_ID\tAFR\tSAS\tEAS\tEUR\tAMR\tOCN\tMES\n" > $(DATA)/$@/ancestries.txt
# 	cat $(DATA)/$@/topmed.global.ancestry.txt | perl -pe 's/ /\t/g' >> $(DATA)/$@/ancestries.txt

# bcf:
# 	mkdir -p $(DATA)/$@
# 	$(foreach c,$(shell seq 1 22) X,ln -s /net/topmed2/working/gt-release/exchange-area/freeze.9b/minDP0/freeze.9b.chr$(c).pass_and_fail.gtonly.minDP0.bcf $(DATA)/$@/chr$(c).bcf$(NL)) 
# 	$(foreach c,$(shell seq 1 22) X,ln -s /net/topmed2/working/gt-release/exchange-area/freeze.9b/minDP0/freeze.9b.chr$(c).pass_and_fail.gtonly.minDP0.bcf.csi $(DATA)/$@/chr$(c).bcf.csi$(NL)) 

# chrom-sizes:
# 	mkdir -p $(DATA)/$@ && cd $(DATA)/$@ && scp wolverine:/lab/data/reference/human/hg38/topmed/index/STAR/STAR_genome_GRCh38_noALT_noHLA_noDecoy_ERCC_v30_oh100/chrNameLength.txt hg38.chrom_sizes

# genotype-concordance:
# 	mkdir -p $(DATA)/$@
# 	zcat /net/topmed2/working/jweinstk/identity_checking_rnaseq/analysis/rnaseq_genotype_concordances_best_matches_2020_07_28.tsv.gz > $(DATA)/$@/best-matches.txt

unrelated:
	mkdir -p $(DATA)/$@
	cp /net/topmed10/working/porchard/rnaseq/work/unrelated-samples/freeze-beta/results/king-2.2.7/unrelated/unrelated.txt $(DATA)/$@/

genotype-pca:
	mkdir -p $(DATA)/$@
	cp /net/topmed10/working/porchard/rnaseq/work/genotype-pca/freeze-beta/unrelated/results/pca/PC-scores.txt $(DATA)/$@/

metadata:
	mkdir -p $(DATA)/$@
	cp /net/topmed10/working/porchard/rnaseq/data/wolverine/metadata.txt $(DATA)/$@/

pass-qc:
	mkdir -p $(DATA)/$@
	cat /net/topmed10/working/porchard/rnaseq/work/qc/all-tissues/results/get-pca-outliers/*.nonoutliers.txt > $(DATA)/$@/pass-qc.txt

gene-expression:
	mkdir -p $(DATA)/$@
	ln -s /net/topmed10/working/porchard/rnaseq/data/wolverine/all-cohorts.gene_counts.hdf5 $(DATA)/$@/counts.hdf5
	ln -s /net/topmed10/working/porchard/rnaseq/data//wolverine/all-cohorts.gene_tpm.hdf5 $(DATA)/$@/tpm.hdf5

mappability:
	mkdir -p $(DATA)/$@ && cd $(DATA)/$@ && wget https://ndownloader.figshare.com/files/13514759 -O hg38_gene_mappability.txt.gz

remapped:
	ln -s /net/topmed10/working/porchard/rnaseq/work/remap-for-sqtl-pass-only $(DATA)/$@

##### ANALYSES #####
eqtl-scan-in-lung: ANALYSIS=$(WORK)/eqtl-scan-in/freeze-beta/Lung
eqtl-scan-in-lung:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && python $(BIN)/select-samples-for-scan.py --unrelated $(DATA)/unrelated/unrelated.txt --in-genotype-pca $(DATA)/genotype-pca/PC-scores.txt --metadata $(DATA)/metadata/metadata.txt --pass-pca-qc $(DATA)/pass-qc/pass-qc.txt --tissues Lung --cohorts LTRC > $(ANALYSIS)/samples-for-scan.txt

eqtl-scan-in-PBMC: ANALYSIS=$(WORK)/eqtl-scan-in/freeze-beta/PBMC
eqtl-scan-in-PBMC:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && python $(BIN)/select-samples-for-scan.py --unrelated $(DATA)/unrelated/unrelated.txt --in-genotype-pca $(DATA)/genotype-pca/PC-scores.txt --metadata $(DATA)/metadata/metadata.txt --pass-pca-qc $(DATA)/pass-qc/pass-qc.txt --tissues PBMC --cohorts MESA > $(ANALYSIS)/samples-for-scan.txt

eqtl-scan-in-monocyte: ANALYSIS=$(WORK)/eqtl-scan-in/freeze-beta/Monocyte
eqtl-scan-in-monocyte:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && python $(BIN)/select-samples-for-scan.py --unrelated $(DATA)/unrelated/unrelated.txt --in-genotype-pca $(DATA)/genotype-pca/PC-scores.txt --metadata $(DATA)/metadata/metadata.txt --pass-pca-qc $(DATA)/pass-qc/pass-qc.txt --tissues Monocyte --cohorts MESA > $(ANALYSIS)/samples-for-scan.txt

eqtl-scan-in-T_cell: ANALYSIS=$(WORK)/eqtl-scan-in/freeze-beta/T_cell
eqtl-scan-in-T_cell:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && python $(BIN)/select-samples-for-scan.py --unrelated $(DATA)/unrelated/unrelated.txt --in-genotype-pca $(DATA)/genotype-pca/PC-scores.txt --metadata $(DATA)/metadata/metadata.txt --pass-pca-qc $(DATA)/pass-qc/pass-qc.txt --tissues T_cell --cohorts MESA > $(ANALYSIS)/samples-for-scan.txt

eqtl-scan-in-whole-blood: ANALYSIS=$(WORK)/eqtl-scan-in/freeze-beta/Whole_blood
eqtl-scan-in-whole-blood:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && python $(BIN)/select-samples-for-scan.py --unrelated $(DATA)/unrelated/unrelated.txt --in-genotype-pca $(DATA)/genotype-pca/PC-scores.txt --metadata $(DATA)/metadata/metadata.txt --pass-pca-qc $(DATA)/pass-qc/pass-qc.txt --tissues Whole_blood --cohorts COPDGene FHS_Levy FHS_Ramachandran SPIROMICS GALAII SAGE WHI > $(ANALYSIS)/samples-for-scan.txt

eqtl-scan-in-nasal-epithelial: ANALYSIS=$(WORK)/eqtl-scan-in/freeze-beta/Nasal_epithelial
eqtl-scan-in-nasal-epithelial:
	mkdir -p $(ANALYSIS)
	cd $(ANALYSIS) && python $(BIN)/select-samples-for-scan.py --unrelated $(DATA)/unrelated/unrelated.txt --in-genotype-pca $(DATA)/genotype-pca/PC-scores.txt --metadata $(DATA)/metadata/metadata.txt --pass-pca-qc $(DATA)/pass-qc/pass-qc.txt --tissues Nasal_epithelial --cohorts COPDGene > $(ANALYSIS)/samples-for-scan.txt

tensorqtl-in-eQTL: ANALYSIS=$(WORK)/tensorqtl-in/eqtl
tensorqtl-in-eQTL:
	mkdir -p $(ANALYSIS)/data
	zcat $(DATA)/gtf/gencode.v30.GRCh38.ERCC.genes.collapsed_only.gtf.gz > $(ANALYSIS)/data/topmed.gtf
	$(foreach t,Lung Nasal_epithelial Whole_blood T_cell PBMC Monocyte,cp $(WORK)/eqtl-scan-in/freeze-beta/$(t)/samples-for-scan.txt $(ANALYSIS)/data/$(t).samples.txt$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --mappability $(DATA)/mappability/hg38_gene_mappability.txt.gz --samples_glob '$(ANALYSIS)/data/*.samples.txt' --metadata $(DATA)/metadata/metadata.txt --tpm $(DATA)/gene-expression/tpm.hdf5 --counts $(DATA)/gene-expression/counts.hdf5 --genotype_pca $(DATA)/genotype-pca/PC-scores.txt --gtf $(ANALYSIS)/data/topmed.gtf $(ROOT)/eqtl-input.nf &

long-indels:
	mkdir -p $(ANALYSIS)
	python $(BIN)/get-long-indels.py $(WORK)/subset-topmed-bcf/freeze-alpha/results/plink-pass-filter/*.bim > $(ANALYSIS)/long-indels.txt

cluster-introns-no-nunique-mask:
	mkdir -p $(ANALYSIS)/data/samples
	mkdir -p $(ANALYSIS)/data/counts
	printf "torid\tnwdid\n" > $(ANALYSIS)/data/tor-to-nwd.txt
	cut-name tor,wgs $(DATA)/metadata/metadata.txt | grep NWD >> $(ANALYSIS)/data/tor-to-nwd.txt # TODO: cut-name
	cp $(WORK)/tensorqtl-in/freeze-beta/data/*.samples.txt $(ANALYSIS)/data/samples/
	zcat $(DATA)/gtf/gencode.v30.GRCh38.ERCC.genes.collapsed_only.gtf.gz > $(ANALYSIS)/data/gencode.v30.GRCh38.ERCC.genes.collapsed_only.gtf
	cd $(ANALYSIS) && nohup nextflow run -resume -qs 1000 --results $(ANALYSIS)/results --leafcutter_dir $(ROOT)/bin --tor_to_nwd $(ANALYSIS)/data/tor-to-nwd.txt --samples_glob '$(ANALYSIS)/data/samples/*' --junction_files_glob '$(DATA)/remapped/results/exon-exon-junction-counts-filtered/*' --collapsed_gtf $(ANALYSIS)/data/gencode.v30.GRCh38.ERCC.genes.collapsed_only.gtf $(ROOT)/cluster-introns-no-nunique-mask.nf &

tensorqtl-in-sQTL: ANALYSIS=$(WORK)/tensorqtl-in/sqtl
tensorqtl-in-sQTL:
	mkdir -p $(ANALYSIS)/data
	$(foreach t,Lung Nasal_epithelial Whole_blood T_cell PBMC Monocyte,cp $(WORK)/cluster-introns-no-nunique-mask/data/samples/$(t).samples.txt $(ANALYSIS)/data/$(t).samples.txt$(NL))
	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --samples_glob '$(ANALYSIS)/data/*.samples.txt' --phenotypes_glob '$(WORK)/cluster-introns-no-nunique-mask/results/ratios/*.leafcutter.bed.gz' --metadata $(DATA)/metadata/metadata.txt --genotype_pca $(DATA)/genotype-pca/PC-scores.txt $(ROOT)/sqtl-input.nf &


# # downsampling saturation analyses
# saturation-samples: ANALYSIS=$(WORK)/saturation/samples
# saturation-samples:
# 	mkdir -p $(ANALYSIS)
# 	cd $(ANALYSIS) && python $(BIN)/make-nested-subsets.py --input $(WORK)/freezes/freeze-1RNA/files/sample-ref/samples/Whole_blood.samples.txt --n $(shell seq 500 500 6500) --prefix Whole_blood.

# saturation-cluster-introns: ANALYSIS=$(WORK)/saturation/cluster-introns
# saturation-cluster-introns:
# 	#mkdir -p $(ANALYSIS)/data/samples
# 	#mkdir -p $(ANALYSIS)/data/counts
# 	#printf "torid\tnwdid\n" > $(ANALYSIS)/data/tor-to-nwd.txt
# 	#cut-name tor,wgs $(DATA)/wolverine/metadata.txt | grep NWD >> $(ANALYSIS)/data/tor-to-nwd.txt
# 	$(foreach t,$(shell seq 500 500 6500),#cp $(WORK)/saturation/samples/Whole_blood.$(t).txt $(ANALYSIS)/data/samples/Whole_blood_$(t).samples.txt$(NL))
# 	#zcat $(DATA)/gtf/gencode.v30.GRCh38.ERCC.genes.collapsed_only.gtf.gz > $(ANALYSIS)/data/gencode.v30.GRCh38.ERCC.genes.collapsed_only.gtf
# 	cd $(ANALYSIS) && nohup nextflow run -resume -qs 1000 -work-dir /net/topmed11/working/porchard/rnaseq/work/$@/work --results $(ANALYSIS)/results --tor_to_nwd $(ANALYSIS)/data/tor-to-nwd.txt --samples_glob '$(ANALYSIS)/data/samples/*' --junction_files_glob '$(WORK)/remap-for-sqtl-pass-only/results/exon-exon-junction-counts-filtered/*' --collapsed_gtf $(ANALYSIS)/data/gencode.v30.GRCh38.ERCC.genes.collapsed_only.gtf $(ROOT)/cluster-introns-no-nunique-mask.nf &

# saturation-tensorqtl-in-eqtl: ANALYSIS=$(WORK)/saturation/tensorqtl-in/eqtl
# saturation-tensorqtl-in-eqtl:
# 	mkdir -p $(ANALYSIS)/data
# 	zcat $(DATA)/gtf/gencode.v30.GRCh38.ERCC.genes.collapsed_only.gtf.gz > $(ANALYSIS)/data/topmed.gtf
# 	$(foreach t,$(shell seq 500 500 6500),cp $(WORK)/saturation/samples/Whole_blood.$(t).txt $(ANALYSIS)/data/Whole_blood_$(t).samples.txt$(NL))
# 	cp $(WORK)/genotype-pca/freeze-beta/unrelated/results/pca/PC-scores.txt $(ANALYSIS)/data/genotype-pcs.txt
# 	cp $(DATA)/wolverine/metadata.txt $(ANALYSIS)/data/metadata.txt
# 	ln -s $(DATA)/wolverine/all-cohorts.gene_counts.hdf5 $(ANALYSIS)/data/counts.hdf5
# 	ln -s $(DATA)/wolverine/all-cohorts.gene_tpm.hdf5 $(ANALYSIS)/data/tpm.hdf5
# 	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --mappability $(DATA)/mappability/hg38_gene_mappability.txt.gz --samples_glob '$(ANALYSIS)/data/*.samples.txt' --metadata $(ANALYSIS)/data/metadata.txt --tpm $(ANALYSIS)/data/tpm.hdf5 --counts $(ANALYSIS)/data/counts.hdf5 --genotype_pca $(ANALYSIS)/data/genotype-pcs.txt --gtf $(ANALYSIS)/data/topmed.gtf $(ROOT)/eqtl-scan-input.nf &

# saturation-tensorqtl-in-sqtl: ANALYSIS=$(WORK)/saturation/tensorqtl-in/sqtl
# saturation-tensorqtl-in-sqtl:
# 	mkdir -p $(ANALYSIS)/data
# 	cp $(WORK)/saturation/cluster-introns/data/samples/*.samples.txt $(ANALYSIS)/data/
# 	cp $(WORK)/genotype-pca/freeze-beta/unrelated/results/pca/PC-scores.txt $(ANALYSIS)/data/genotype-pcs.txt
# 	cp $(DATA)/wolverine/metadata.txt $(ANALYSIS)/data/metadata.txt
# 	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --samples_glob '$(ANALYSIS)/data/*.samples.txt' --phenotypes_glob '$(WORK)/saturation/cluster-introns/results/ratios/*.leafcutter.bed.gz' --metadata $(ANALYSIS)/data/metadata.txt --genotype_pca $(ANALYSIS)/data/genotype-pcs.txt $(ROOT)/tensorqtl-in-sqtl.nf &

# # ancestry-specific e/sQTL input
# ancestry-specific-samples: ANALYSIS=$(WORK)/ancestry-specific/samples
# ancestry-specific-samples: TISSUES=Lung Nasal_epithelial T_cell PBMC Monocyte Whole_blood
# ancestry-specific-samples:
# 	mkdir -p $(ANALYSIS)
# 	$(foreach t,$(TISSUES),python $(BIN)/assign-sample-ancestry.py --fraction 0.75 /net/topmed10/working/porchard/rnaseq/work/freezes/freeze-1RNA/files/sample-ref/samples/$(t).samples.txt /net/topmed10/working/porchard/rnaseq/data/wolverine/metadata.txt /net/topmed10/working/porchard/rnaseq/data/ancestry/ancestries.txt > $(ANALYSIS)/$(t).samples.txt$(NL))
# 	mkdir -p $(ANALYSIS)/split-by-ancestry
# 	$(foreach t,$(TISSUES),grep EUR $(ANALYSIS)/$(t).samples.txt | cut -f1 > $(ANALYSIS)/split-by-ancestry/$(t)___EUR.samples.txt$(NL))
# 	$(foreach t,PBMC Whole_blood,grep AFR $(ANALYSIS)/$(t).samples.txt | cut -f1 > $(ANALYSIS)/split-by-ancestry/$(t)___AFR.samples.txt$(NL))
# 	$(foreach t,PBMC,grep EAS $(ANALYSIS)/$(t).samples.txt | cut -f1 > $(ANALYSIS)/split-by-ancestry/$(t)___EAS.samples.txt$(NL))

# ancestry-specific-cluster-introns: ANALYSIS=$(WORK)/ancestry-specific/cluster-introns
# ancestry-specific-cluster-introns:
# 	mkdir -p $(ANALYSIS)/data/samples
# 	mkdir -p $(ANALYSIS)/data/counts
# 	printf "torid\tnwdid\n" > $(ANALYSIS)/data/tor-to-nwd.txt
# 	cut-name tor,wgs $(DATA)/wolverine/metadata.txt | grep NWD >> $(ANALYSIS)/data/tor-to-nwd.txt
# 	cp $(WORK)/ancestry-specific/samples/split-by-ancestry/*.samples.txt $(ANALYSIS)/data/samples/
# 	zcat $(DATA)/gtf/gencode.v30.GRCh38.ERCC.genes.collapsed_only.gtf.gz > $(ANALYSIS)/data/gencode.v30.GRCh38.ERCC.genes.collapsed_only.gtf
# 	cd $(ANALYSIS) && nohup nextflow run -resume -qs 1000 -work-dir /net/topmed3/working/porchard/rnaseq/work/$@/work --results $(ANALYSIS)/results --tor_to_nwd $(ANALYSIS)/data/tor-to-nwd.txt --samples_glob '$(ANALYSIS)/data/samples/*' --junction_files_glob '$(WORK)/remap-for-sqtl-pass-only/results/exon-exon-junction-counts-filtered/*' --collapsed_gtf $(ANALYSIS)/data/gencode.v30.GRCh38.ERCC.genes.collapsed_only.gtf $(ROOT)/cluster-introns-no-nunique-mask.nf &

# ancestry-specific-tensorqtl-in-eqtl: ANALYSIS=$(WORK)/ancestry-specific/tensorqtl-in/eqtl
# ancestry-specific-tensorqtl-in-eqtl:
# 	mkdir -p $(ANALYSIS)/data
# 	zcat $(DATA)/gtf/gencode.v30.GRCh38.ERCC.genes.collapsed_only.gtf.gz > $(ANALYSIS)/data/topmed.gtf
# 	cp $(WORK)/ancestry-specific/samples/split-by-ancestry/*.samples.txt $(ANALYSIS)/data/
# 	cp $(WORK)/genotype-pca/freeze-beta/unrelated/results/pca/PC-scores.txt $(ANALYSIS)/data/genotype-pcs.txt
# 	cp $(DATA)/wolverine/metadata.txt $(ANALYSIS)/data/metadata.txt
# 	ln -s $(DATA)/wolverine/all-cohorts.gene_counts.hdf5 $(ANALYSIS)/data/counts.hdf5
# 	ln -s $(DATA)/wolverine/all-cohorts.gene_tpm.hdf5 $(ANALYSIS)/data/tpm.hdf5
# 	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --plink_glob '$(WORK)/subset-topmed-bcf/freeze-beta/results/plink/chr*' --mappability $(DATA)/mappability/hg38_gene_mappability.txt.gz --samples_glob '$(ANALYSIS)/data/*.samples.txt' --metadata $(ANALYSIS)/data/metadata.txt --tpm $(ANALYSIS)/data/tpm.hdf5 --counts $(ANALYSIS)/data/counts.hdf5 --genotype_pca $(ANALYSIS)/data/genotype-pcs.txt --gtf $(ANALYSIS)/data/topmed.gtf $(ROOT)/scan-input-freeze-alpha.nf &

# ancestry-specific-tensorqtl-in-sqtl: ANALYSIS=$(WORK)/ancestry-specific/tensorqtl-in/sqtl
# ancestry-specific-tensorqtl-in-sqtl:
# 	mkdir -p $(ANALYSIS)/data
# 	cp $(WORK)/ancestry-specific/cluster-introns/data/samples/*.samples.txt $(ANALYSIS)/data/
# 	cp $(WORK)/genotype-pca/freeze-beta/unrelated/results/pca/PC-scores.txt $(ANALYSIS)/data/genotype-pcs.txt
# 	cp $(DATA)/wolverine/metadata.txt $(ANALYSIS)/data/metadata.txt
# 	cd $(ANALYSIS) && nohup nextflow run -resume --results $(ANALYSIS)/results --samples_glob '$(ANALYSIS)/data/*.samples.txt' --phenotypes_glob '$(WORK)/ancestry-specific/cluster-introns/results/ratios/*.leafcutter.bed.gz' --metadata $(ANALYSIS)/data/metadata.txt --genotype_pca $(ANALYSIS)/data/genotype-pcs.txt $(ROOT)/tensorqtl-in-sqtl.nf &
