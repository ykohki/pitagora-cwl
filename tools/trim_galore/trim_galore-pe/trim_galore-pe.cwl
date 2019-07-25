class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  sbg: 'https://www.sevenbridges.com'
id: trim_galore
baseCommand:
  - trim_galore
inputs:
  - id: read1
    type: File
    inputBinding:
      position: 0
  - id: read2
    type: File
    inputBinding:
      position: 0
  - id: fastqc
    type: boolean?
    inputBinding:
      position: -2
      prefix: '--fastqc'
  - id: trim1
    type: boolean?
    inputBinding:
      position: -1
      prefix: '--trim1'
  - id: paired
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--paired'
outputs:
  - id: out1
    type: File
    outputBinding:
      glob: $(inputs.read1+'_1_val_1.fq.gz')
  - id: out2
    type: File
    outputBinding:
      glob: $(inputs.read1+'_2_val_2.fq.gz')
label: trim_galore
requirements:
  - class: DockerRequirement
    dockerPull: miasteinberg/trim-galore
  - class: InlineJavascriptRequirement
