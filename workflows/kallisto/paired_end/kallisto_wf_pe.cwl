class: Workflow
cwlVersion: v1.0
$namespaces:
  edam: 'http://edamontology.org/'
  s: 'https://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: bootstrap_samples
    type: int?
    'sbg:x': 598.8079833984375
    'sbg:y': 326.895263671875
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
    'sbg:x': -95.8778076171875
    'sbg:y': 44.586036682128906
  - id: nthreads
    type: int?
    'sbg:x': 771.1320190429688
    'sbg:y': 502.8089599609375
  - id: index_name
    type: string
    'sbg:x': 323.92694091796875
    'sbg:y': -359.7203369140625
  - id: paired
    type: boolean?
    'sbg:x': 113.95063781738281
    'sbg:y': 256.9856872558594
  - id: fastqc
    type: boolean?
    'sbg:x': 338.44940185546875
    'sbg:y': 352.3373107910156
  - id: trim1
    type: boolean?
    'sbg:x': 118.77857208251953
    'sbg:y': 430.7911682128906
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
    'sbg:x': 282.7830505371094
    'sbg:y': -101.06982421875
  - id: forward
    outputSource:
      - fasterq_dump/forward
    type: File?
    'sbg:x': 410.37158203125
    'sbg:y': -38.82793045043945
  - id: fastqFiles
    outputSource:
      - fasterq_dump/fastqFiles
    type: 'File[]'
    'sbg:x': 374.0997619628906
    'sbg:y': 64.37157440185547
  - id: index_file
    outputSource:
      - kallisto_index/index_file
    type: File
    'sbg:x': 678.4219970703125
    'sbg:y': -270.0132751464844
  - id: out2
    outputSource:
      - trim_galore/out2
    type: File
    'sbg:x': 835.7261962890625
    'sbg:y': -179.94200134277344
  - id: out1
    outputSource:
      - trim_galore/out1
    type: File
    'sbg:x': 944.3546142578125
    'sbg:y': 439.24005126953125
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
      - id: nthreads
        source: nthreads
      - id: out_dir_name
        source: out_dir_name
    out:
      - id: quant_output
    run: >-
      https://raw.githubusercontent.com/pitagora-network/pitagora-cwl/master/tools/kallisto/quant/paired_end/kallisto_quant_pe.cwl
    label: 'kallisto quant: runs the quantification algorithm'
    'sbg:x': 969.80859375
    'sbg:y': 82.40431213378906
  - id: kallisto_index
    in:
      - id: fasta_files
        source:
          - fasta_files
      - id: index_name
        source: index_name
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
    'sbg:x': 160.92019653320312
    'sbg:y': 38.55112075805664
  - id: trim_galore
    in:
      - id: read1
        source: fasterq_dump/forward
      - id: read2
        source: fasterq_dump/reverse
      - id: fastqc
        source: fastqc
      - id: trim1
        source: trim1
      - id: paired
        source: paired
    out:
      - id: out1
      - id: out2
    run: ../../../tools/trim_galore/trim_galore-pe/trim_galore-pe.cwl
    label: trim_galore
    'sbg:x': 559.8977661132812
    'sbg:y': 143.93017578125
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
