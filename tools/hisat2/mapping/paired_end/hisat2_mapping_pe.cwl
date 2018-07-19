cwlVersion: v1.0
class: CommandLineTool
label: "HISAT2: graph-based alignment of next generation sequencing reads to a population of genomes [cwl for paired end read]"
doc: "HISAT2: graph-based alignment of next generation sequencing reads to a population of genomes. https://ccb.jhu.edu/software/hisat2/manual.shtml#command-line"

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/hisat2:2.1.0--py36h2d50403_1

baseCommand: [hisat2]

arguments:
  - prefix: -S
    valueFrom: $(runtime.outdir)/$(inputs.out_sam_name)
  - prefix: -x
    valueFrom: $(inputs.hisat2_idx_basedir.path)/$(inputs.hisat2_idx_basename)

inputs:
  hisat2_idx_basedir:
    label: "Path to the directory the index for the reference genome are in"
    doc: "Path to the directory the index for the index files, such as .1.ht2 / etc exist."
    type: Directory
  hisat2_idx_basename:
    label: "Basename of the hisat2 index files"
    doc: "Basename of the hisat2 index files, not including extensions like .1.ht2"
    type: string
  fq1:
    label: "Read file containing mate 1s"
    doc: "Comma-separated list of files containing mate 1s (filename usually includes _1), e.g. -1 flyA_1.fq,flyB_1.fq. Sequences specified with this option must correspond file-for-file and read-for-read with those specified in <m2>. Reads may be a mix of different lengths."
    type: File
    inputBinding:
      prefix: "-1"
  fq2:
    label: "Read file containing mate 2s"
    doc: "Comma-separated list of files containing mate 2s (filename usually includes _2), e.g. -2 flyA_2.fq,flyB_2.fq. Sequences specified with this option must correspond file-for-file and read-for-read with those specified in <m1>. Reads may be a mix of different lengths. "
    type: File
    inputBinding:
      prefix: "-2"
  out_sam_name:
    label: "Name of file to write SAM alignments to"
    doc: "Name of file to write SAM alignments to (default: out.sam)"
    type: string?
    default: out.sam
  dta:
    label: "Report alignments tailored for transcript assemblers"
    doc: "Report alignments tailored for transcript assemblers including StringTie. With this option, HISAT2 requires longer anchor lengths for de novo discovery of splice sites. This leads to fewer alignments with short-anchors, which helps transcript assemblers improve significantly in computation and memory usage."
    type: boolean?
    default: true
    inputBinding:
      prefix: --downstream-transcriptome-assembly
  dta_cufflinks:
    label: "Report alignments tailored specifically for Cufflinks"
    doc: "Report alignments tailored specifically for Cufflinks. In addition to what HISAT2 does with the above option (--dta), With this option, HISAT2 looks for novel splice sites with three signals (GT/AG, GC/AG, AT/AC), but all user-provided splice sites are used irrespective of their signals. HISAT2 produces an optional field, XS:A:[+-], for every spliced alignment"
    type: boolean?
    default: true
    inputBinding:
      prefix: --dta-cufflinks
  time:
    label: "Print the wall-clock time"
    doc: "Print the wall-clock time required to load the index files and align the reads. This is printed to the 'standard error' ('stderr') filehandle."
    type: boolean?
    default: true
    inputBinding:
      prefix: --time
  nthreads:
    label: "Launch NTHREADS parallel search threads"
    doc: "Launch NTHREADS parallel search threads (default: 1). Threads will run on separate processors/cores and synchronize when parsing reads and outputting alignments. Searching for alignments is highly parallel, and speedup is close to linear. Increasing -p increases HISAT2's memory footprint. E.g. when aligning to a human genome index, increasing -p from 1 to 8 increases the memory footprint by a few hundred megabytes. This option is only available if bowtie is linked with the pthreads library (i.e. if BOWTIE_PTHREADS=0 is not specified at build time)."
    type: int
    inputBinding:
      prefix: --threads

outputs:
  hisat2_sam:
    type: File
    outputBinding:
      glob: $(inputs.out_sam_name)
