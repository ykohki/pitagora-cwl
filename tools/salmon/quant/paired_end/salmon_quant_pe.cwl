class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  edam: 'http://edamontology.org/'
  s: 'https://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
baseCommand:
  - salmon
  - quant
inputs:
  - id: fq1
    type: File
    inputBinding:
      position: 0
      prefix: '--mates1'
    label: 'File containing the #1 mates'
    doc: 'File containing the #1 mates'
  - id: fq2
    type: File
    inputBinding:
      position: 0
      prefix: '--mates2'
    label: 'File containing the #2 mates'
    doc: 'File containing the #2 mates'
  - id: gc_bias
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--gcBias'
    label: '[beta for single-end reads] Perform fragment GC bias correction'
    doc: '[beta for single-end reads] Perform fragment GC bias correction'
  - id: gene_map
    type: File?
    inputBinding:
      position: 0
      prefix: '--geneMap'
    label: File containing a mapping of transcripts to genes
    doc: >-
      File containing a mapping of transcripts to genes. If this file is
      provided salmon will output both quant.sf and quant.genes.sf files, where
      the latter contains aggregated gene-level abundance estimates. The
      transcript to gene mapping should be provided as either a GTF file, or a
      in a simple tab-delimited format where each line contains the name of a
      transcript and the gene to which it belongs separated by a tab. The
      extension of the file is used to determine how the file should be parsed.
      Files ending in '.gtf', '.gff' or '.gff3' are assumed to be in GTF format;
      files with any other extension are assumed to be in the simple format. In
      GTF / GFF format, the 'transcript_id' is assumed to contain the transcript
      identifier and the 'gene_id' is assumed to contain the corresponding gene
      identifier.
  - default: 0
    id: gibbs_samples
    type: int
    inputBinding:
      position: 0
      prefix: '--numGibbsSamples'
    label: Number of Gibbs sampling rounds to perform.
    doc: Number of Gibbs sampling rounds to perform.
  - default: 0
    id: incompat_prior
    type: int
    inputBinding:
      position: 0
      prefix: '--incompatPrior'
    label: >-
      the prior probability that an alignment that disagrees with the specified
      library type (--libType) results from the true fragment origin
    doc: >-
      This option sets the prior probability that an alignment that disagrees
      with the specified library type (--libType) results from the true fragment
      origin. Setting this to 0 specifies that alignments that disagree with the
      library type should be 'impossible', while setting it to 1 says that
      alignments that disagree with the library type are no less likely than
      those that do
  - id: index_dir
    type: Directory
    inputBinding:
      position: 0
      prefix: '--index'
    label: salmon index
    doc: salmon index
  - default: IU
    id: lib_type
    type: string
    inputBinding:
      position: 0
      prefix: '--libType'
    label: Format string describing the library type
    doc: Format string describing the library type
  - default: false
    id: meta
    type: boolean
    inputBinding:
      position: 0
      prefix: '--meta'
    label: a metagenomic dataset
    doc: >-
      If you're using Salmon on a metagenomic dataset, consider setting this
      flag to disable parts of the abundance estimation model that make less
      sense for metagenomic data.
  - default: 2
    id: nthreads
    type: int
    inputBinding:
      position: 0
      prefix: '--threads'
    label: The number of threads to use concurrently.
    doc: The number of threads to use concurrently.
  - default: 0
    id: num_bootstraps
    type: int
    inputBinding:
      position: 0
      prefix: '--numBootstraps'
    label: Number of bootstrap samples to generate
    doc: >-
      Number of bootstrap samples to generate. Note: This is mutually exclusive
      with Gibbs sampling.
  - default: salmon_quant
    id: quant_out_dirname
    type: string
    label: Output quantification directory
    doc: Output quantification directory.
  - id: seq_bias
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--seqBias'
    label: Perform sequence-specific bias correction.
    doc: Perform sequence-specific bias correction.
  - id: validateMappings
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--validateMappings'
outputs:
  - id: quant_results
    type: Directory
    outputBinding:
      glob: $(inputs.quant_out_dirname)
doc: >-
  Salmon is a tool for quantifying the expression of transcripts using RNA-seq
  data. Salmon uses new algorithms (specifically, coupling the concept of
  quasi-mapping with a two-phase inference procedure) to provide accurate
  expression estimates very quickly (i.e. wicked-fast) and while using little
  memory. Salmon performs its inference using an expressive and realistic model
  of RNA-seq data that takes into account experimental attributes and biases
  commonly observed in real RNA-seq data.
  http://salmon.readthedocs.io/en/latest/ (no documentation for command line
  options, see salmon quant --help-reads. salmon quant has many advanced
  options, here are only basic options defined.)
label: 'Salmon quant: quantifying the samples'
arguments:
  - position: 0
    prefix: '--output'
    valueFrom: $(runtime.outdir)/$(inputs.quant_out_dirname)
hints:
  - class: DockerRequirement
    dockerPull: 'combinelab/salmon:0.10.2'
requirements:
  - class: InlineJavascriptRequirement
$schemas:
  - 'https://schema.org/docs/schema_org_rdfa.html'
  - 'http://edamontology.org/EDAM_1.18.owl'
's:author':
  - class: 's:Person'
    's:email': 'mailto:inutano@gmail.com'
    's:identifier': 'https://orcid.org/0000-0003-3777-5945'
    's:name': Tazro Ohta
's:codeRepository': 'https://github.com/pitagora-network/pitagora-cwl'
's:license': 'https://spdx.org/licenses/Apache-2.0'
