####Activate Qiime2 environment 
conda activate qiime2-2021.8

####Prepare the manifest file. The manifest file is a tab-delimited file that links the samples name to the fastq files (full paths).


qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path manifest_file.tsv --output-path paired-end-demux.qza --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize --i-data paired-end-demux.qza --o-visualization demux.qzv

qiime dada2 denoise-paired --i-demultiplexed-seqs paired-end-demux.qza --p-trunc-len-f 240 --p-trunc-len-r 53 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza

qiime phylogeny align-to-tree-mafft-fasttree --i-sequences rep-seqs.qza --o-alignment aligned-rep-seqs.qza --o-masked-alignment masked-aligned-rep-seqs.qza --o-tree unrooted-tree.qza --o-rooted-tree rooted-tree.qza

####Download the pre-trained clasifier Silva 138 99% OTUs from 515F/806R region of sequences at https://docs.qiime2.org/2022.2/data-resources/
qiime feature-classifier classify-sklearn --i-classifier silva-138-99-515-806-nb-classifier.qza --i-reads rep-seqs.qza --o-classification taxonomy.qza

mkdir phyloseq
qiime tools export --input-path rep-seqs.qza --output-path phyloseq

qiime tools export --input-path unrooted-tree.qza --output-path phyloseq
cd phyloseq
mv tree.nwk unrooted_tree.nwk
cd ../

qiime tools export --input-path rooted-tree.qza --output-path phyloseq
cd phyloseq
mv tree.nwk rooted_tree.nwk
cd ../

qiime metadata tabulate --m-input-file taxonomy.qza --o-visualization taxonomy.qzv

qiime tools export --input-path table.qza --output-path phyloseq

biom convert —-i feature-table.biom -o phyloseq/otu_table.txt --to-tsv

qiime tools export --input-path taxonomy.qza --output-path phyloseq/

