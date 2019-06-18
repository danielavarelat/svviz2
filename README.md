# svviz2

This is a forked repo from (https://github.com/nspies/svviz2). 

Modified to get just the read count information in a report per sample as svviz2 does not natively support parallelization.

Installation
------------
svviz2 requires **python 3.3** or greater. 


Documentation
-------------

More in-depth documentation of the original svviz2 is available at [https://svviz2.readthedocs.io](https://svviz2.readthedocs.io).

Usage
-----

```
ssw library not found
usage: svviz2 [options] --ref REF --variants VARIANTS BAM [BAM2 ...]

svviz2 version 2.0a3

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  bam                   sorted, indexed bam file containing reads of interest to plot; can be specified multiple
                        times to load multiple samples
  --ref REF, -r REF     reference fasta file (a .faidx index file will be created if it doesn't exist so you need
                        write permissions for this directory)
  --variants VARIANTS, -V VARIANTS
                        the variants to analyze, in vcf or bcf format (vcf files may be compressed with gzip)

Optional arguments:
  --outdir OUTDIR, -o OUTDIR
                        output directory for visualizations, summaries, etc (default: current working directory)
  --savereads           output the read realignments against the appropriate alt or ref allele (default: false)
  --min-mapq MIN_MAPQ   only reads with mapq>=MIN_MAPQ will be analyzed; when analyzing paired-end data,
                        at least one read end must be near the breakpoints with this mapq (default:0)
  --align-distance ALIGN_DISTANCE
                        sequence upstream and downstream of breakpoints to include when performing re-alignment
                        (default: infer from data)
  --batch-size BATCH_SIZE
                        Number of reads to analyze at once; larger batch-size values may run more quickly
                        but will require more memory (default=10000)
  --downsample DOWNSAMPLE
                        Ensure the total number of reads per event per sample does not exceed this number
                        by downsampling (default: infinity)
  --aligner ALIGNER     The aligner to use for realigning reads; either ssw (smith-waterman) or
                        bwa (default=bwa)
  --only-realign-locally
                        Only when using bwa as the aligner backend, when this option is enabled,
                        reads will only be aligned locally around the breakpoints and not also against
                        the full reference genome (default: False)
  --fast                More aggressively skip reads that are unlikely to overlap
                        the breakpoints (default: false)
  --first-variant FIRST_VARIANT
                        Skip all variants before this variant; counting starts with first variant
                        in input VCF as 0 (default: 0)
  --last-variant LAST_VARIANT
                        Skip all variants after this variant; counting starts with first variant
                        in input VCF as 0 (default: end of vcf)

  --report
  
  --no-report
```


```
Example command for just one sample:

    svviz2 {~/sample.bam} 
    --ref {~/gr37.fasta}
    --variants {~/file.vcf}
    --outdir {outdir} 
    --report
```

Recommendation
-----
If the installation present problems with installation, it is recommended to run this modified tool using a built container for https://github.com/papaemmelab/toil_circosigv.

```
Command with singularity:
    singularity exec --bind /home,/ifs,/work  {singularity_image.img} svviz2 ...
```
 
