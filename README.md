# CCBR ChIP-Seq Pipeline


## Installation

1. install nextflow
	```bash
	curl -fsSL get.nextflow.io | bash
	```
	For details, visit https://www.nextflow.io/

2. clone this git repository
	```bash
	#for latest experimental features
	git clone https://github.com/CCBR/ChIP-Seq-Pipeline.git
	```
	
	Otherwise,
	download from the following link
	https://github.com/CCBR/ChIP-Seq-Pipeline/releases/latest
	

## How to run the pipeline

```bash
./nextflow run ChIP-Seq-Pipeline/main.nf \
	-config ChIP-Seq-Pipeline/config \
	--macsconfig='examples/macs.config' \
	--reads='ChIP-Seq-Pipeline/examples/*.fastq' \
	--genome='hg19'
```
	
* Note that the reads names are supposed to end with *".fastq"* or *".fastq.gz"*.
If you experience an error message related to missing file(s), this might be the reason.

* Config file is designed to work at our NIH biowulf2 HPC. If you plan to use on helix, add -profile 'local' option at the end of command line arguments. If you plan to use outside NIH, you need to edit the "genome" section according to your paths information.

* macs.config is a config file for MACS2 to find peaks. ChIP sample and control (input or IgG) with their label concatenated by comma form a line, e.g. `chip_sample_id,control_id,sample_name`. `chip_sample_id` or `control_id` are the basename of the FASTQ files wihtout the trailing *".fastq"* or *".fastq.gz"*.

* Currently, the pipeline fully supports hg19, mm10 and partially supports GRCh38 and GRCm38. Contact us for other genomes. 


#### Additional options: 
* `-resume`                         --> to resume the previous failed run
* `-profile 'local'`                --> for running tools "locally" not thourgh high performance computer queueing mechanisms.
* `-with-timeline 'timeline.html'`  --> record the run time of each process.
* `-with-dag 'flowchart.png'`       --> draw the flowchart
* `--outdir 'ChIP-Seq-Pipeline-Output'`       --> specify the directory to save directory from the default output

We implemented the pipeline using Nextflow.


## Implemented tools
Thanks to the authors of the tools!

1. Trimgalor
2. BWA mem
3. Picard (MarkDuplicate)
4. FASTQC
5. DeepTools
6. Macs2
7. Sicer
8. MEME-ChIP
9. PhantomPeakqualTool
10. CEAS

 
 
### Todo list

1. ChipSeeker
2. Homer
3. PeakDiff
4. ngsplot
 
 
 
### Thanks

I got many nice implementation ideas from Nextflow examples, especially from NGI-ChIP-seq pipeline.
Many thanks to the NGI pipeline developers and Alexei, Ashley, George, Parthav and other CCBR team mebers.
 
 
#### Questions or Suggestions

Email to Bong-Hyun.Kim at NIH dot GOV.

