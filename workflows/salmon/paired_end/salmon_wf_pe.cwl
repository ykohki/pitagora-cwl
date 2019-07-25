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
  - id: libType
    type: string?
    'sbg:x': 416.8983154296875
    'sbg:y': 275.1154479980469
  - id: validateMappings
    type: boolean?
    'sbg:x': 1020.0255126953125
    'sbg:y': -117.04666137695312
  - id: quant_out_dirname
    type: string
    'sbg:x': 850.4768676757812
    'sbg:y': -181.75987243652344
  - id: nthreads
    type: int
    'sbg:x': 265.4693908691406
    'sbg:y': 131.45208740234375
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
  - id: quant_results
    outputSource:
      - salmon_quant_pe/quant_results
    type: Directory
    'sbg:x': 1316.4119873046875
    'sbg:y': 157.33737182617188
  - id: salmon_index_1
    outputSource:
      - salmon_index/salmon_index
    type: Directory
    'sbg:x': 854.864501953125
    'sbg:y': 566.3367919921875
  - id: out2
    outputSource:
      - trim_galore/out2
    type: File
    'sbg:x': 779.7973022460938
    'sbg:y': -172.68836975097656
  - id: out1
    outputSource:
      - trim_galore/out1
    type: File
    'sbg:x': 750.0292358398438
    'sbg:y': 137.93507385253906
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
  - id: salmon_quant_pe
    in:
      - id: fq1
        source: trim_galore/out1
      - id: fq2
        source: trim_galore/out2
      - id: index_dir
        source: salmon_index/salmon_index
      - id: nthreads
        source: nthreads
      - id: quant_out_dirname
        source: quant_out_dirname
      - id: libType
        source: libType
      - id: validateMappings
        source: validateMappings
    out:
      - id: quant_results
    run: ../../../tools/salmon/quant/paired_end/salmon_quant_pe.cwl
    label: 'Salmon quant: quantifying the samples'
    'sbg:x': 1000.611328125
    'sbg:y': 156.0433349609375
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
    'sbg:x': 724.143798828125
    'sbg:y': 0.7432571649551392
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
