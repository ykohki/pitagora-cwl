class: Workflow
cwlVersion: v1.0
$namespaces:
  edam: 'http://edamontology.org/'
  s: 'https://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: bootstrap_samples
    type: int?
    'sbg:x': 390.2039489746094
    'sbg:y': 345.4114990234375
  - id: out_dir_name
    type: string?
    'sbg:x': 557.6082763671875
    'sbg:y': 432.0071716308594
  - id: fasta_files
    type: 'File[]'
    'sbg:x': 272.02874755859375
    'sbg:y': -236.0035858154297
  - id: srafile
    type: File
    'sbg:x': 125.33853912353516
    'sbg:y': 47.872798919677734
outputs:
  - id: quant_output
    outputSource:
      - kallisto_quant/quant_output
    type: Directory
    'sbg:x': 1329.3809814453125
    'sbg:y': 87.38993835449219
  - id: reverse
    outputSource:
      - fasterq_dump/reverse
    type: File?
    'sbg:x': 465.04022216796875
    'sbg:y': -68.56204986572266
  - id: forward
    outputSource:
      - fasterq_dump/forward
    type: File?
    'sbg:x': 569.4714965820312
    'sbg:y': -67.3616943359375
  - id: fastqFiles
    outputSource:
      - fasterq_dump/fastqFiles
    type: 'File[]'
    'sbg:x': 426.3917236328125
    'sbg:y': 173.80682373046875
steps:
  - id: kallisto_quant
    in:
      - id: bootstrap_samples
        source: bootstrap_samples
      - id: fq1
        source: trim_galore/out1
      - id: fq2
        source: trim_galore/out2
      - id: index_file
        source: kallisto_index/index_file
      - id: out_dir_name
        source: out_dir_name
    out:
      - id: quant_output
    run: >-
      https://raw.githubusercontent.com/pitagora-network/pitagora-cwl/master/tools/kallisto/quant/paired_end/kallisto_quant_pe.cwl
    label: 'kallisto quant: runs the quantification algorithm'
    'sbg:x': 941.5993041992188
    'sbg:y': 68.80682373046875
  - id: kallisto_index
    in:
      - id: fasta_files
        source:
          - fasta_files
    out:
      - id: index_file
    run: ../../../tools/kallisto/index/kallisto_index.cwl
    label: >-
      kallisto index: builds an index from a FASTA formatted file of target
      sequences
    'sbg:x': 516.4833374023438
    'sbg:y': -214.7996368408203
  - id: fasterq_dump
    in:
      - id: srafile
        source: srafile
    out:
      - id: fastqFiles
      - id: forward
      - id: reverse
    run: ../../../tools/fasterq-dump/fasterq-dump.cwl
    label: >-
      fasterq-dump: dump .sra format file to generate fastq file, way more
      faster
    'sbg:x': 311.19317626953125
    'sbg:y': 57
  - id: trim_galore
    in:
      - id: read1
        source: fasterq_dump/forward
      - id: read2
        source: fasterq_dump/reverse
    out:
      - id: out1
      - id: out2
    run: ../../../tools/trim_galore/trim_galore-pe/trim_galore-pe.cwl
    label: trim_galore
    'sbg:x': 658.2980346679688
    'sbg:y': 109.09112548828125
requirements: []
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
