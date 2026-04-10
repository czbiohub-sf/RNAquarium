# RNAquarium pipeline parameters



## Reference genome options

Reference genome related files and options required for the workflow.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `genomeSize` |  | `number` |  | True |  |
| `starRefIndexes` |  | `string` |  |  |  |
| `starRefIndexesErcc` |  | `string` |  |  |  |
| `hisatRefIndexes` |  | `string` |  |  |  |
| `refGenome` |  | `string` | Danio_rerio.GRCz11.dna_sm.primary_assembly.fa |  |  |
| `refGenomeGtf` |  | `string` | Danio_rerio.GRCz11.108.gtf |  |  |
| `erccFa` |  | `string` | ERCC92.fa |  |  |
| `erccGtf` |  | `string` | ERCC92.gtf |  |  |

## Input/output options

Basic inputs and outputs

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `accessionList` |  | `string` |  | True |  |
| `fastqPath` |  | `string` |  |  |  |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `help` | Display help text. | `boolean` |  |  | True |
| `h` | Display help text. | `boolean` |  |  |  |
| `reportdir` | directory to write nextflow reports/traces to | `string` | reports |  |  |
| `timestamp` | timestamp to append to nextflow reports (default: current time, YYMMDD-HHmmss) | `string` | 231108-222816 |  | True |
| `tmp` |  | `string` |  |  |  |
| `parallelDownloads` |  | `number` |  |  |  |
| `publishDir` |  | `string` |  |  |  |
| `publishIntermediate` |  | `boolean` |  |  |  |
| `publishFastqs` |  | `boolean` | True |  |  |
| `publishPricefiltered` |  | `boolean` | True |  |  |
| `publishReadcounts` |  | `boolean` | True |  |  |
| `publishHisat` |  | `boolean` | True |  |  |
| `publishStar` |  | `boolean` | True |  |  |

## Pipeline options

Control for pipeline flow

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `skipHostCounts` |  | `boolean` |  |  |  |
| `skipHisat` |  | `boolean` |  |  |  |
| `hisatUseTranscript` |  | `boolean` |  |  |  |

## Tool extra parameters



| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `starIndexGenOptions` |  | `string` |  |  |  |
| `hisatIndexGenOptions` |  | `string` |  |  |  |
| `sraPrefetchOptions` |  | `string` |  |  |  |
| `fastqDumpOptions` |  | `string` |  |  |  |
| `fastpOptions` |  | `string` |  |  |  |
| `priceOptions` |  | `string` |  |  |  |
| `starCountOptions` |  | `string` |  |  |  |
| `samtoolsSortOptions` |  | `string` |  |  |  |
| `htseqCountOptions` |  | `string` |  |  |  |
| `hisatOptions` |  | `string` |  |  |  |
| `starOptions` |  | `string` |  |  |  |
