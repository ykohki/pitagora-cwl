cwlVersion: v1.0
class: CommandLineTool

hints:
  DockerRequirement:
    dockerPull: quay.io/inutano/sra-toolkit:v2.9.0

baseCommand: [pfastq-dump]

inputs:
  sraFiles:
    type: File[]
    inputBinding:
      position: 50
  nthreads:
    type: int
    inputBinding:
      prefix: -t
  split_files:
    type: boolean?
    default: true
    inputBinding:
      prefix: --split-files
  split_spot:
    type: boolean?
    default: true
    inputBinding:
      prefix: --split-spot
  skip_technical:
    type: boolean?
    default: true
    inputBinding:
      prefix: --skip-technical
  readids:
    type: boolean?
    default: true
    inputBinding:
      prefix: --readids
  gzip:
    type: boolean?
    default: true
    inputBinding:
      prefix: --gzip

outputs:
  fastqFiles:
    type: File[]
    outputBinding:
      glob: "*fastq*"
  forward:
    type: File?
    outputBinding:
      glob: "*_1.fastq*"
  reverse:
    type: File?
    outputBinding:
      glob: "*_2.fastq*"