class: Workflow
cwlVersion: v1.0
$namespaces:
  edam: 'http://edamontology.org/'
  s: 'https://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: transcript_fasta
    type: File
    'sbg:x': 256.9574279785156
    'sbg:y': 514.9254760742188
  - id: srafile
    type: File
    'sbg:x': 138.68988037109375
    'sbg:y': -115.03089141845703
  - id: nthreads_1
    type: int
    'sbg:x': 469.1220703125
    'sbg:y': 702.6729125976562
  - id: index_name
    type: string
    'sbg:x': 445.0498962402344
    'sbg:y': 577.9353637695312
  - id: paired
    type: boolean?
    'sbg:x': 111.96560668945312
    'sbg:y': 275.1271057128906
  - id: fastqc
    type: boolean?
    'sbg:x': 251.74615478515625
    'sbg:y': 373.4911804199219
  - id: trim1
    type: boolean?
    'sbg:x': 434.2374267578125
    'sbg:y': -314.05609130859375
  - id: validateMappings
    type: boolean?
    'sbg:x': 1014.980712890625
    'sbg:y': -122.70207977294922
  - id: lib_type
    type: string
    'sbg:x': 1163.4071044921875
    'sbg:y': 574.6221313476562
  - id: nthreads
    type: int
    'sbg:x': 1005.178955078125
    'sbg:y': 696.4437866210938
  - id: quant_out_dirname
    type: string
    'sbg:x': 883.3572387695312
    'sbg:y': -226.32052612304688
outputs:
  - id: reverse
    outputSource:
      - fasterq_dump/reverse
    type: File?
    'sbg:x': 559.8676147460938
    'sbg:y': -216.40318298339844
  - id: forward
    outputSource:
      - fasterq_dump/forward
    type: File?
    'sbg:x': 703.4783325195312
    'sbg:y': -229.67811584472656
  - id: fastqFiles
    outputSource:
      - fasterq_dump/fastqFiles
    type: 'File[]'
    'sbg:x': 530.904052734375
    'sbg:y': 49.095664978027344
  - id: salmon_index_1
    outputSource:
      - salmon_index/salmon_index
    type: Directory
    'sbg:x': 854.864501953125
    'sbg:y': 566.3367919921875
  - id: out1
    outputSource:
      - trim_galore/out1
    type: File
    'sbg:x': 876.8758544921875
    'sbg:y': 329.1932678222656
  - id: out2
    outputSource:
      - trim_galore/out2
    type: File
    'sbg:x': 792.7486572265625
    'sbg:y': -53.908966064453125
  - id: quant_results
    outputSource:
      - salmon_quant_pe/quant_results
    type: Directory
    'sbg:x': 1394.4482421875
    'sbg:y': 202.15577697753906
steps:
  - id: salmon_index
    in:
      - id: index_name
        source: index_name
      - id: nthreads
        source: nthreads_1
      - id: transcript_fasta
        source: transcript_fasta
    out:
      - id: salmon_index
    run: ../../../tools/salmon/index/salmon_index.cwl
    label: 'Salmon index: building an index'
    'sbg:x': 681.2225952148438
    'sbg:y': 436.9718933105469
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
    'sbg:x': 425.3819885253906
    'sbg:y': -108.58637237548828
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
    'sbg:x': 724.942626953125
    'sbg:y': 84.00498962402344
  - id: salmon_quant_pe
    in:
      - id: fq1
        source: trim_galore/out1
      - id: fq2
        source: trim_galore/out2
      - id: index_dir
        source: salmon_index/salmon_index
      - id: lib_type
        source: lib_type
      - id: nthreads
        source: nthreads
      - id: quant_out_dirname
        source: quant_out_dirname
      - id: validateMappings
        source: validateMappings
    out:
      - id: quant_results
    run: ../../../tools/salmon/quant/paired_end/salmon_quant_pe.cwl
    label: 'Salmon quant: quantifying the samples'
    'sbg:x': 1090.2020263671875
    'sbg:y': 202.59974670410156
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
