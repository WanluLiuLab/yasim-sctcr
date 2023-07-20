all: sim_tcr_50.fq sim_notcr_50.d

.SECONDARY: hg38.ncbiRefSeq.gtf hg38.ncbiRefSeq_chr7_14.gtf hg38.fai hg38.fa

hg38.ncbiRefSeq.gtf.gz:
	axel https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz

hg38.ncbiRefSeq.gtf: hg38.ncbiRefSeq.gtf.gz
	gzip -dcf hg38.ncbiRefSeq.gtf.gz > hg38.ncbiRefSeq.gtf

hg38.ncbiRefSeq_chr7_14.gtf: hg38.ncbiRefSeq.gtf
	grep -e '^chr7\s' -e '^chr14\s' hg38.ncbiRefSeq.gtf > hg38.ncbiRefSeq_chr7_14.gtf

hg38.fa.gz:
	axel https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz

hg38.fa: hg38.fa.gz
	gzip -dcf hg38.fa.gz > hg38.fa

hg38.fai: hg38.fa
	samtools faidx hg38.fa

distclean: clean
	rm -rf hg38*

clean:
	rm -fr tcr* sim* notcr*

tcr_barcode.txt:
	python -m yasim_sc generate_barcode -n 10 -o tcr_barcode.txt

tcr_depth.tsv: tcr_barcode.txt
	python -m yasim_sctcr generate_tcr_depth \
		-b tcr_barcode.txt \
		-o tcr_depth.tsv \
		-d 20

tcr_cache.json: hg38.ncbiRefSeq_chr7_14.gtf hg38.fai
	python -m yasim_sctcr generate_tcr_cache \
		--tcr_genelist_path ../data/tcr_genelist.min.json.xz \
		-o tcr_cache.json \
		--tcr_aa_table_path ../data/IMGT_Protein_Display.min.json.xz \
		-f hg38.fa \
		-g hg38.ncbiRefSeq_chr7_14.gtf

sim_tcr.nt.fa: tcr_cache.json tcr_barcode.txt
	python -m yasim_sctcr rearrange_tcr \
		--tcr_genelist_path ../data/tcr_genelist.min.json.xz \
		--tcr_cache_path tcr_cache.json \
		--cdr3_deletion_table_path ../data/cdr3_deletion_table.min.json.xz \
		--cdr3_insertion_table_path ../data/cdr3_insertion_table.min.json.xz \
		-b tcr_barcode.txt \
		-o sim_tcr
	python -m labw_utils.bioutils split_fasta sim_tcr.nt.fa

sim_tcr_50.fq: tcr_depth.tsv sim_tcr.nt.fa
	python -m yasim art \
		-F sim_tcr.nt.fa.d \
		-o sim_tcr_50 \
		--sequencer_name GA2 \
		--read_length 50 \
		-d tcr_depth.tsv \
		-e art_illumina \
		-j 20

hg38.ncbiRefSeq_subsampled.gtf: hg38.ncbiRefSeq_chr7_14.gtf
	python -m yasim_sc sample_pcg \
		-i ../data/ncbi_dataset.tsv.xz \
		-g hg38.ncbiRefSeq_chr7_14.gtf \
		-o hg38.ncbiRefSeq_subsampled.gtf \
		--num_genes_to_sample 20

notcr_gene_depth.tsv: hg38.ncbiRefSeq_subsampled.gtf
	python -m yasim generate_gene_depth \
		-g hg38.ncbiRefSeq_subsampled.gtf \
		-d 5 \
		-o notcr_gene_depth.tsv

notcr_isoform_depth.tsv: notcr_gene_depth.tsv hg38.ncbiRefSeq_subsampled.gtf
	python -m yasim generate_isoform_depth \
		-g hg38.ncbiRefSeq_subsampled.gtf \
		-d notcr_gene_depth.tsv \
		-o notcr_isoform_depth.tsv

notcr_isoform_depth_sc.d: notcr_isoform_depth.tsv
	python -m yasim_sc generate_barcoded_isoform_replicates \
		-d notcr_isoform_depth.tsv \
		-b tcr_barcode.txt \
		-o notcr_isoform_depth_sc.d

notcr_trans.fa: hg38.ncbiRefSeq_subsampled.gtf hg38.fai
	python -m labw_utils.bioutils transcribe \
		-f hg38.fa \
		-g hg38.ncbiRefSeq_subsampled.gtf \
		-o notcr_trans.fa

sim_notcr_50.d: notcr_trans.fa notcr_isoform_depth_sc.d
	python -m yasim_sc art \
		-F notcr_trans.fa.d \
		-o sim_notcr_50.d \
		-d notcr_isoform_depth_sc.d \
		-e art_illumina \
		-j 20
