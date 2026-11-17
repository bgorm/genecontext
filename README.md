# genecontext

[![Build and test](https://github.com/bgorm/genecontext/actions/workflows/build.yml/badge.svg)](https://github.com/bgorm/genecontext/actions/workflows/build.yml)

Annotate VCF/BCF variants with nearby genes using a GTF or GFF3 annotation.
**Genecontext is a streaming [bcftools](https://github.com/samtools/bcftools) plugin:** it adds INFO fields while
preserving the original variants, samples, and FORMAT fields.
It works with general VCF/BCF files and is designed particularly for annotating GWAS loci in [GWAS-VCF](https://github.com/MRCIEU/gwas-vcf-specification)
workflows using bcftools plugins
from the [score suite](https://github.com/freeseek/score) and [gwaslocusfinder](https://github.com/bgorm/gwaslocusfinder).

Choose nearest genes and TSSs, distances, gene neighborhoods, promoter
overlaps, and exon/CDS/intron annotations.

- [Usage](#usage)
- [Installation](#installation)
- [Example: GWAS summary statistics to annotated loci](#example-gwas-summary-statistics-to-annotated-loci)
- [Annotation examples](#annotation-examples)
- [Gene annotation files](#gene-annotation-files)
- [Annotation definitions](#annotation-definitions)
- [INFO fields](#info-fields)
- [Performance and limitations](#performance-and-limitations)
- [License](#license)

## Usage

```sh
bcftools +genecontext input.bcf -Ob -o annotated.bcf -W -- \
    -a genes.gtf.gz
```

The defaults add the nearest gene, nearest gene-boundary TSS, and their
distances: `GC_NEAREST`, `GC_DISTANCE`, `GC_TSS`, and `GC_TSS_DISTANCE`.
All gene biotypes are included. Nearest searches stop at 1 Mb from POS.

General bcftools options go before `--`; genecontext options go after it.
Use `-Ou` to stream uncompressed BCF into another bcftools command.
The complete plugin help is shown below; the [INFO table](#info-fields)
describes the output groups.

<!-- BEGIN GENECONTEXT HELP -->
```text
About: Annotate variant gene neighborhoods from GTF/GFF3 in a bcftools stream.
(genecontext version 0.2.2)

Usage: bcftools +genecontext [General Options] <in.vcf|in.bcf|-> -- [Plugin Options]

   bcftools handles -i/-e, -r/-R, -t/-T, -o, -O, --threads and -W before --.
   All new annotations are INFO fields and use the variant's 1-based POS.

Plugin options (after --):
   -a, --annotation FILE      Gene annotation: plain/gzip/BGZF GTF or GFF3 (required)
       --emit LIST            Comma-separated output groups [nearest,distance,nearest-tss,tss-distance]
                              nearest, distance, nearest-tss, tss-distance, neighbors,
                              neighbors-tss, map, overlap, ids, biotype, strand,
                              flanks, count, overlap-count, tss-details, promoter,
                              exon, cds, intron, feature-details, relations; or all
       --tss-source gene|canonical  Gene boundaries or tagged transcript TSSs [gene]
       --canonical-tag STR    Transcript tag used for canonical selection [Ensembl_canonical]
       --promoter-upstream INT   Promoter bp upstream of the selected TSS [2000]
       --promoter-downstream INT Promoter bp downstream of the selected TSS [200]
       --feature-transcripts all|canonical  Transcripts used for exon/CDS/intron output [all]
       --max-distance INT     Maximum nearest-gene/TSS/flank distance in bp [1000000]
       --window INT           Neighborhood radius on each side of POS in bp [250000]
       --dot-bp INT           Base pairs per dot in the map, rounded up [10000]
       --max-dots INT         Maximum dots per map gap; ~ marks longer gaps [80]
       --label name|id        Gene labels; missing names always fall back to IDs [name]
       --protein-coding       Shortcut for --biotypes protein_coding
       --biotypes LIST        Include exact comma-separated gene biotypes [all]
       --exclude-biotypes LIST  Exclude exact comma-separated gene biotypes
       --features LIST        Gene-level feature types [gene,pseudogene,ncRNA_gene]
       --transcript-features LIST  Transcript feature types [transcript,mRNA and common RNA types]
       --annotation-format auto|gtf|gff3  Attribute syntax [auto]
       --chr-map FILE         Two columns: VCF contig, annotation contig [exact names]
       --missing-contig MODE  warn, error or ignore [warn]
       --prefix STR           Prefix for emitted INFO tags [GC_]
       --overwrite            Replace selected existing tags of compatible Number/Type
       --sort-mem INT         Approximate sort buffer in MiB for preparation [64]
       --temp-dir DIR         Scratch directory parent [TMPDIR or system temporary directory]
       --quiet                Suppress preparation and completion messages
   -h, --help                 Show this help

Definitions:
   Gene body includes both endpoints. Distance is zero inside the body. All ties survive.
   TSS uses the gene or tagged transcript start on + and end on -. Unknown strands have no TSS.
   Canonical mode has no fallback. All transcripts bearing the chosen tag are eligible.
   Promoter offsets are strand-aware; introns are gaps between merged exons of one transcript.
   Transcript features are required for canonical TSSs; exon/CDS records supply feature overlaps.
   Neighborhood distances are negative to the genomic left, positive to the right.
   Map anchors are the nearest body boundary or POS for [overlapping|genes]; [] is intergenic.
   Dots measure gaps between anchors and window edges. Tied anchors use |.
   Existing compatible CSI/TBI indexes are reused; otherwise a temporary table is prepared.
   Missing nearest matches remain missing; searches do not extend beyond max-distance.

Examples:
   bcftools +genecontext index.bcf -Ob -o genes.bcf -W -- -a genes.gtf.gz
   bcftools view -Ou input.bcf | bcftools +genecontext -Oz -o genes.vcf.gz -- \
       -a genes.gff3.gz --emit all --protein-coding --window 100000 --dot-bp 10000
   bcftools +genecontext index.bcf -- -a full.gtf.gz --tss-source canonical \
       --emit nearest-tss,tss-distance,tss-details,promoter,exon,cds,intron,feature-details
```
<!-- END GENECONTEXT HELP -->

## Installation

Requirements: bcftools, a C++17 compiler, GNU make, and matching HTSlib
headers.

### Default build

Clone the repository and build the plugin:

```sh
git clone https://github.com/bgorm/genecontext.git
cd genecontext
bash scripts/bootstrap.sh
export BCFTOOLS_PLUGINS="$PWD/build:${BCFTOOLS_PLUGINS:-}"
bcftools +genecontext --help
```

The bootstrap downloads HTSlib 1.25 headers if needed and builds
`build/genecontext.so`. Use bcftools/HTSlib 1.25 with this default build.

### Other bcftools installations

To build against another bcftools/HTSlib installation, supply matching headers:

```sh
make -j4 BCFTOOLS=/path/to/bcftools PLUGIN_HTSLIB=/path/to/htslib-source
```

For installed headers whose version cannot be obtained with `pkg-config`:

```sh
make -j4 PLUGIN_HTSLIB=/path/to/include PLUGIN_HTSLIB_VERSION=1.25
```

Run `make clean` before changing the host version, compiler, or build flags.
Standard `CPPFLAGS`, `CXXFLAGS`, `LDFLAGS`, and `LDLIBS` are supported. The
plugin checks bcftools and HTSlib versions when loaded and resolves HTSlib
symbols from its host. Builds and integration tests with bcftools/HTSlib 1.25
have passed on Linux and macOS in [GitHub Actions](https://github.com/bgorm/genecontext/actions/runs/44137618439).

### Install and test

An explicit path works without setting `BCFTOOLS_PLUGINS`:
`bcftools +/path/to/genecontext.so ...`.
To install the plugin in a persistent directory:

```sh
make install plugindir="$HOME/.local/libexec/bcftools"
export BCFTOOLS_PLUGINS="$HOME/.local/libexec/bcftools:${BCFTOOLS_PLUGINS:-}"
```

`DESTDIR` is supported for staging. Run `make test` for the synthetic
integration suite; tests require Python 3, bcftools, bgzip, and tabix.

The [build workflow](.github/workflows/build.yml) runs on Linux and macOS with
bcftools/HTSlib 1.25. It tests the plugin, checks installation, and repeats the
integration suite in a fresh public source export.

### Worked example dependencies

The worked example also needs `bcftools +munge` from
[freeseek/score](https://github.com/freeseek/score#installation),
[gwaslocusfinder](https://github.com/bgorm/gwaslocusfinder#installation),
and the commands `wget`, `gzip`, `awk`, GNU `sort`, `bgzip`, and `tabix`.
Make both additional plugins available on the bcftools plugin search path.

## Example: GWAS summary statistics to annotated loci

This example uses the **GRCh38 Alzheimer's disease summary statistics** from
[Bellenguez et al. (2022)](https://doi.org/10.1038/s41588-022-01024-z), GWAS Catalog accession
[GCST90027158](https://www.ebi.ac.uk/gwas/studies/GCST90027158).
The specific file is also used in the
[score conversion examples](https://github.com/freeseek/score/blob/master/examples.md#alzheimer).
All reference files and coordinates below are GRCh38.

Run the following steps in order in a Bash shell. Downloads include a full
reference genome and genome-wide summary statistics; allow several GB of
disk space.

### 1. Download the summary statistics and reference genome

```sh
mkdir -p genecontext-example
cd genecontext-example

wget -O GCST90027158_buildGRCh38.tsv.gz \
    https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/GCST90027158_buildGRCh38.tsv.gz

(
    set -euo pipefail
    wget -qO- https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz |
        gzip -dc > GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
)
```

This is the **GRCh38 no-ALT analysis set** used in
[Giulio Genovese's score manual](https://github.com/freeseek/score#installation),
following [Heng Li's recommendation](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use).
It uses `chr` chromosome names, matching the GENCODE annotation below.
The FASTA lets `+munge` identify the reference allele and orient effect
estimates to ALT. `+munge` creates its `.fai` index when needed.

### 2. Download, sort, and index the gene annotation

Use the [GENCODE v50 comprehensive chromosome annotation](https://www.gencodegenes.org/human/release_50.html):

```sh
(
    set -euo pipefail
    wget -qO- https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gtf.gz |
        gzip -dc |
        awk '!/^#/ && NF' |
        LC_ALL=C sort -S 256M -t $'\t' -k1,1 -k4,4n -k5,5n |
        bgzip -@ 4 -c > gencode.v50.sorted.gtf.gz
    tabix -f -C -p gff gencode.v50.sorted.gtf.gz
)
```

This saves a reusable GTF and its `.csi` index. It retains transcript, exon,
and CDS records so the same file supports every annotation option.

### 3. Convert the summary statistics to indexed GWAS-VCF

```sh
(
    set -euo pipefail
    bcftools +munge -c SSF -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
        -s AD_2022 -Ou \
        GCST90027158_buildGRCh38.tsv.gz |
        bcftools view -f PASS,. -m2 -M2 -Ou |
        bcftools sort -Ob -o ad.bcf -W
)
```

This writes `ad.bcf` and its index. The `SSF` preset maps the summary-statistic
columns, including the P-value to `FORMAT/LP` (−log10 P). The filter keeps
biallelic records that passed the reference-allele checks.
Gwaslocusfinder requires this sorted, indexed input file.

### 4. Find significant loci and stream them into genecontext

```sh
(
    set -euo pipefail
    bcftools +gwaslocusfinder -Ou ad.bcf -- \
        --p1 5e-8 --distance-kb 1000 |
        bcftools +genecontext -Ob -o ad.loci.genes.bcf -W -- \
            -a gencode.v50.sorted.gtf.gz \
            --emit nearest,distance --max-distance 1000000 \
            --missing-contig error
)
```

Gwaslocusfinder selects variants with P ≤ 5 × 10⁻⁸, joins consecutive
significant positions separated by at most 1 Mb, and keeps the strongest
variant from each locus. These **locus index variants** pass directly into
genecontext through uncompressed BCF (`-Ou`). For LD-clumped lead variants,
use the [LD reference options](https://github.com/bgorm/gwaslocusfinder#reference-panels)
in gwaslocusfinder.

The output, `ad.loci.genes.bcf`, contains the index variants, their locus
annotations, and two new INFO fields:

| Field | Meaning |
| --- | --- |
| `GC_NEAREST` | Nearest gene body; all equally close genes are retained |
| `GC_DISTANCE` | Distance in bp; zero when POS is inside the gene |

This annotates each index variant's position. The locus boundaries remain
available in `FORMAT/LF_START` and `FORMAT/LF_END`.

### 5. Export the results

```sh
bcftools query -H \
    -f '%CHROM\t%POS\t%ID[\t%LP\t%LF_START\t%LF_END]\t%INFO/GC_NEAREST\t%INFO/GC_DISTANCE\n' \
    ad.loci.genes.bcf > ad.loci.genes.tsv

head ad.loci.genes.tsv
```

The table includes chromosome, index position, variant ID, −log10 P, locus
start/end, nearest gene, and distance. A missing nearest match is `.`.

## Annotation examples

These examples reuse the locus variants and GTF from the worked example.
`--emit` selects the output groups; `--emit all` selects every group.

### Protein-coding genes

Add every annotation, considering only protein-coding genes:

```sh
bcftools +genecontext ad.loci.genes.bcf -Ob -o ad.loci.coding.bcf -W -- \
    -a gencode.v50.sorted.gtf.gz --protein-coding --emit all --overwrite
```

`--overwrite` refreshes the nearest-gene fields already present in the input.
For another selection, replace `--protein-coding` with, for example,
`--biotypes protein_coding,lncRNA` or `--exclude-biotypes processed_pseudogene`.

### Gene neighborhoods and maps

List genes within 250 kb and draw a neighborhood map with 10 kb per dot:

```sh
bcftools +genecontext ad.loci.genes.bcf -Ob -o ad.loci.neighborhood.bcf -W -- \
    -a gencode.v50.sorted.gtf.gz \
    --emit neighbors,neighbors-tss,map --window 250000 --dot-bp 10000

bcftools query -H -f '%CHROM\t%POS\t%INFO/GC_NEIGHBORS\t%INFO/GC_MAP\n' \
    ad.loci.neighborhood.bcf > ad.loci.neighborhood.tsv
```

Neighborhood distances are negative to the genomic left and positive to the
right. In a map such as `..GENEA..[GENEB]...GENEC..`, brackets identify the
gene containing POS; `[]` marks an intergenic variant. Dots show approximate
gaps between genes' nearest boundaries and the variant.

### Canonical TSSs, promoters, and transcript overlaps

Use transcripts tagged `Ensembl_canonical` for TSSs and promoters, and add
exon/CDS/intron annotations:

```sh
bcftools +genecontext ad.loci.genes.bcf -Ob -o ad.loci.transcripts.bcf -W -- \
    -a gencode.v50.sorted.gtf.gz --tss-source canonical \
    --emit nearest-tss,tss-distance,tss-details,promoter,exon,cds,intron,feature-details \
    --promoter-upstream 2000 --promoter-downstream 200
```

Exon/CDS/intron annotations use all transcripts by default. Add
`--feature-transcripts canonical` to restrict them to the same tag, or
`--canonical-tag MANE_Select` to select MANE transcripts. Genes without the
chosen tag have no canonical TSS; there is no gene-boundary fallback.

### Select samples and filter variants

Select the `AD_2022` sample and records whose FILTER is `PASS` or missing,
then stream into the plugin. With a multi-sample input, replace `AD_2022`
with the phenotype to retain:

```sh
(
    set -euo pipefail
    bcftools view -s AD_2022 -f PASS,. -Ou ad.loci.genes.bcf |
        bcftools +genecontext -Oz -o ad.loci.filtered.context.vcf.gz -W -- \
            -a gencode.v50.sorted.gtf.gz \
            --emit nearest,distance,neighbors,map --overwrite
)
```

### Keep separate annotation sets

Use different prefixes to retain annotations for all genes and protein-coding
genes in the same output. This example emits the default four fields for
each selection, under `ALL_` and `PC_`:

```sh
(
    set -euo pipefail
    bcftools +genecontext ad.loci.genes.bcf -Ou -- \
        -a gencode.v50.sorted.gtf.gz --prefix ALL_ |
        bcftools +genecontext -Ob -o ad.loci.two-sets.bcf -W -- \
            -a gencode.v50.sorted.gtf.gz \
            --protein-coding --prefix PC_
)
```

### MANE transcript annotations

```sh
bcftools +genecontext ad.loci.genes.bcf -Ob -o ad.loci.mane.bcf -W -- \
    -a gencode.v50.sorted.gtf.gz --protein-coding \
    --tss-source canonical --canonical-tag MANE_Select \
    --feature-transcripts canonical \
    --emit nearest-tss,tss-distance,tss-details,promoter,exon,cds,intron,feature-details
```

## Gene annotation files

Plain or compressed GTF/GFF3 files are accepted. An existing tabix index is
reused; otherwise the plugin prepares a temporary index automatically.

For the equivalent GENCODE GFF3 download:

```sh
(
    set -euo pipefail
    {
        printf '##gff-version 3\n'
        wget -qO- https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gff3.gz |
            gzip -dc |
            awk '/^##FASTA([[:space:]]|$)/ { fasta=1 }
                 !fasta && !/^#/ && NF' |
            LC_ALL=C sort -S 256M -t $'\t' -k1,1 -k4,4n -k5,5n
    } | bgzip -@ 4 -c > gencode.v50.sorted.gff3.gz
    tabix -f -C -p gff gencode.v50.sorted.gff3.gz
)

bcftools +genecontext ad.loci.genes.bcf -Ob -o ad.loci.gff3.bcf -W -- \
    -a gencode.v50.sorted.gff3.gz --emit nearest,distance --overwrite
```

Use matching genome builds and chromosome names. GENCODE's **CHR** files
cover the standard chromosomes; **PRI** adds primary-assembly scaffolds;
**ALL** also includes patches and alternative loci. PRI and CHR give the
same annotations on their shared chromosomes for the same release and
annotation set. See [chromosome mapping](#chromosome-names)
if your VCF uses names such as `1` instead of `chr1`.

### Select gene biotypes

Using the locus variants from the worked example, include protein-coding
and long noncoding RNA genes:

```sh
bcftools +genecontext ad.loci.genes.bcf -Ob -o ad.loci.coding-lncRNA.bcf -W -- \
    -a gencode.v50.sorted.gtf.gz \
    --biotypes protein_coding,lncRNA --overwrite
```

Or consider all gene biotypes except processed pseudogenes:

```sh
bcftools +genecontext ad.loci.genes.bcf -Ob -o ad.loci.filtered-genes.bcf -W -- \
    -a gencode.v50.sorted.gtf.gz \
    --exclude-biotypes processed_pseudogene --overwrite
```

Biotype matching is exact and case-sensitive. `--protein-coding` is shorthand
for `--biotypes protein_coding`; supply only one of these. Exclusions can be
combined with inclusions. Gene records with no biotype are eligible when all
types are selected, but cannot match an explicit inclusion list. No biotype
is inferred from CDS/exon records.

### Feature types and attributes

The default feature types are `gene,pseudogene,ncRNA_gene`. Override them with
`--features` for providers using different gene-level features. The plugin
does not reconstruct genes from transcript/exon-only annotations.

Default transcript feature types include:

- `transcript`, `mRNA`, `ncRNA`, `lnc_RNA`, `lncRNA`, `tRNA`, `rRNA`, `miRNA`,
  `pre_miRNA`, `snRNA`, `snoRNA`, `scRNA`, `primary_transcript`.
- `pseudogenic_transcript`, `V_gene_segment`, `D_gene_segment`, `J_gene_segment`,
  `C_gene_segment`, `VD_gene_segment`.
- `RNA`, `lincRNA`, `nc_primary_transcript`, `processed_transcript`,
  `aberrant_processed_transcript`, `NMD_transcript_variant`.

These cover the transcript families in Ensembl's
[GFF3 mapping](https://www.ensembl.info/2017/07/06/gff3-and-sequence-ontology-terms/)
and the transcript-level terms accepted by its
[release 116 GFF reader](https://github.com/Ensembl/ensembl-vep/blob/release/116/modules/Bio/EnsEMBL/VEP/AnnotationSource/File/GFF.pm).
`--transcript-features` **replaces** the entire list; retain any default types
you still need when supplying custom types. Gene and transcript feature types
must be distinct when transcript annotations are requested. `exon` and `CDS`
are reserved for structural features in these modes.

Attributes are normalized as follows, using the first nonempty value:

| Value | Attribute priority |
| --- | --- |
| Gene ID | `gene_id`, `ID`, `locus_tag` |
| Gene name | `gene_name`, `Name`, `gene` |
| Biotype | `gene_type`, `gene_biotype`, `biotype` |
| Transcript ID | `transcript_id`, `ID` on a transcript feature |

GTF attributes and GFF3 `key=value` attributes are detected automatically;
`--annotation-format gtf` or `gff3` forces a format. GTF accepts quoted
strings and unquoted scalar values such as `level 2` and `exon_number 1`. GFF3
percent escapes are decoded before consistent output escaping. The first
value of scalar attributes is used; all `tag` and `Parent` values are retained.
Comments, blank lines, and nonselected feature types are skipped; preparation
stops at a GFF3 `##FASTA` directive. Selected gene/transcript rows with invalid
coordinates, missing identifiers, or malformed attributes produce an error
with source context.

GTF children are associated using `gene_id` and `transcript_id`. GFF3
`Parent` links resolve gene → transcript → exon/CDS relationships, including
multiple parents. Raw GFF3 `ID` aliases are retained when a gene or transcript
also supplies a different display identifier. Percent-encoded commas in IDs
are preserved as part of the ID, not treated as parent-list separators.

Children without a resolvable eligible gene are omitted. GTF exons with
explicit gene/transcript IDs can supply introns even when the file lacks a
transcript row; they cannot supply a tagged transcript TSS or qualify for
`--feature-transcripts canonical`. Transcript-based annotations load the full
annotated spans of candidate genes before resolving their children. Conflicting
or multipart coordinates for the same resolved transcript ID within those
spans produce an error, including when only one part intersects the original
search window. Child features are expected to lie within their parent gene's
annotated span. A missing structural field means no resolved overlap;
this also occurs when the input lacks the necessary child features.

### Chromosome names

Use an annotation release matching the VCF's genome assembly. Contig names
must match exactly; no automatic `chr` stripping or assembly conversion takes
place. An optional two-column `--chr-map` maps VCF contigs to annotation
contigs without renaming the VCF itself. For a sorted input named
`variants.bcf` with unprefixed chromosome names:

```sh
awk 'BEGIN {
    for (i=1; i<=22; i++) printf "%d\tchr%d\n", i, i
    printf "X\tchrX\nY\tchrY\nMT\tchrM\n"
}' > gencode.chr-map.tsv

bcftools +genecontext variants.bcf -Ob -o variants.mapped-genes.bcf -W -- \
    -a gencode.v50.sorted.gtf.gz \
    --chr-map gencode.chr-map.tsv --missing-contig error
```

This maps VCF contigs `1`–`22`, `X`, `Y`, and `MT` to GENCODE's `chr1`–`chr22`,
`chrX`, `chrY`, and `chrM`. Add explicit mappings for other names present in
your input.

Missing contigs warn once by default. `--missing-contig error` makes them
fatal; `--missing-contig ignore` leaves their fields missing silently.

The attribute conventions follow the primary format references from
[GENCODE](https://www.gencodegenes.org/pages/data_format.html) and
[NCBI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/).

## Annotation definitions

### Positions and distances

All calculations use the variant's **POS**, a single 1-based position. They
do not use REF length, INFO/END, allele-specific spans, or phenotype-specific
`FORMAT/LF_START`/`LF_END`. Multiallelic records therefore receive one shared
set of annotations. Input records may be unsorted unless the requested
bcftools operation, such as output indexing, requires sorted input.

For a gene spanning inclusive coordinates `[start,end]`:

- Body distance is `start-POS` to its left, zero within the body, and
  `POS-end` to its right. A variant one base outside a boundary has distance 1.
- The TSS is `start` on `+` and `end` on `−`, using the selected gene or
  transcript feature. Its absolute distance is `abs(TSS-POS)`. Features with
  unknown strand do not supply TSSs or promoters.
- Neighborhood lists use **gene coordinate minus variant coordinate**:
  negative is genomic left, positive is genomic right. For body distances,
  the coordinate is the nearest boundary, or POS when overlapping.
  Signs do not change with gene strand and do not denote transcriptional
  upstream/downstream.

`--max-distance` limits nearest body/TSS/flank searches; `--window` separately
limits both neighborhood lists, the map, relationships, and the neighborhood
count. Limits are inclusive. Zero permits exact overlaps only. Both options
accept integer base pairs up to 2,147,483,647. Fetching complete gene models
does not relax these annotation limits.
Promoter widths have their own limits and do not depend on these two options.
Exon/CDS/intron output tests POS directly, independently of either window.

Tied nearest genes are ordered by annotation chromosome, start, end, and ID,
with remaining attributes breaking ties. Neighborhood lists are ordered by
their signed distance, then gene coordinates. Exact duplicate features are
removed. Multiple rows with the same identifier are treated as separate parts
of one gene: each metric chooses the nearest part and reports the ID once,
without filling any gap between parts. Filters are applied before this step.
Flanks exclude a gene ID from both sides if **any part contains POS**. Otherwise,
the nearest part on each side is eligible; a variant in a gap can have the same
multipart gene as both flanks. TSS and promoter detail fields report each unique
gene ID/TSS position/transcript/strand combination once.

### TSS selection and Open Targets

`--tss-source gene` is the default and uses **gene-feature boundaries**.
`--tss-source canonical` uses the start/end of transcript features
with the exact `Ensembl_canonical` tag. Repeated GTF `tag` attributes and
comma-separated GFF3 `tag` values are both supported. A different tag can
be selected explicitly, for example `--canonical-tag MANE_Select`.
Canonical selection uses the annotation's tags, without ranking transcripts
by length or guessing a substitute. The definitions of these tags are given
by [GENCODE](https://www.gencodegenes.org/pages/tags.html) and
[Ensembl](https://www.ensembl.org/info/genome/genebuild/canonical.html).

Open Targets documents TSS distance relative to the **canonical transcript**.
Its pipeline notes derive that TSS from the transcript's strand and
coordinates. Canonical mode follows this coordinate definition; exact
agreement also requires the same assembly, annotation release, transcript
selection, and gene filters.
[Open Targets L2G documentation](https://platform-docs.opentargets.org/gentropy/locus-to-gene-l2g),
[target/transcript pipeline notes](https://github.com/opentargets/issues/issues/4457).

For example, a plus-strand gene could span 100–1000 while its canonical
transcript starts at 200. At variant POS 150, gene mode uses TSS 100 and
canonical mode uses TSS 200. This is a schematic example, not a claim about
a specific gene.

Genes without a resolvable, tagged transcript with known strand remain
eligible for body annotations but do not supply a canonical TSS. There is
**no fallback to the gene boundary**. A file containing only gene records
therefore cannot provide canonical TSS annotations.

If a custom tag selects several transcripts for one gene, all their TSSs are
eligible. Gene label/ID lists report a gene once at its closest site; the TSS
neighborhood chooses the first site in deterministic gene/position/transcript
order for equal distances. `GC_TSS_DETAIL` retains every globally tied nearest
site, including multiple transcripts of one gene. Its source is `gene` or
`canonical`; the transcript ID is `.` for gene boundaries. Request
`--emit tss-details` when transcript identity matters.

### Promoters and transcript features

`--emit promoter` uses the TSS source selected by `--tss-source`. With
upstream width U and downstream width D, the inclusive promoter interval is
`[TSS-U,TSS+D]` on `+` and `[TSS-D,TSS+U]` on `−`. Defaults are 2,000 bp
upstream and 200 bp downstream. Both widths can be zero and are configured
with `--promoter-upstream` and `--promoter-downstream`.

`GC_PROMOTER_DETAIL` reports the transcript/gene TSS and the variant's
transcription-oriented offset: `(POS-TSS)` on `+`, `(TSS-POS)` on `−`.
Negative means upstream, positive means downstream, and zero is the TSS.
All overlapping selected promoters are retained. These are configurable
positional promoter windows, not experimentally measured promoter activity.

`--emit exon,cds,intron,feature-details` adds gene labels and transcript-level
overlap details. Exon and CDS annotations use the input's explicit feature
intervals. Introns are derived separately for each gene/transcript pair:
overlapping or adjacent exons are merged, then each internal gap is an
intron. Gene/transcript tails outside the annotated exons are not introns.
Single-exon transcripts have no derived introns; CDS intervals alone are not
used to infer them.

A variant can be exonic in one transcript and intronic in another, so the
gene can appear in both `GC_EXON` and `GC_INTRON`. The detailed tuples make
this explicit, for example:

```text
GC_FEATURE=GENEA|G1|T1|intron|1401|1999,GENEA|G1|T2|exon|1350|2050
```

Feature coordinates are 1-based inclusive. `FEATURE` is `exon`, `CDS`, or
`intron`. Direct gene-level exon/CDS records can have `TRANSCRIPT_ID=.`;
these records do not imply introns. Structural outputs are overlap labels,
not variant-effect predictions, and still use POS rather than the full span
of an indel or structural variant.

`--feature-transcripts all` is the default and is independent of
`--tss-source`. To restrict structural annotations to the selected canonical
tag, use `--feature-transcripts canonical`. This also excludes child features
without a resolvable tagged transcript. Biotype filters always use the
**gene** record's biotype, including when an alternative transcript has a
different biotype.

`--emit relations` describes the variant relative to each gene body in the
body neighborhood. Upstream/downstream follow gene strand, unlike the
genomic signs in `GC_NEIGHBORS`. Overlapping variants are `overlap`; unknown
strand outside a body produces `unknown_strand`. These relationships use
gene bodies, independently of TSS selection.

### Neighborhood map

A map might look like:

```text
..GENEA.GENEB..[GENEC]....GENEX..
```

`[GENEC]` means POS falls within GENEC's body. Multiple overlapping genes
appear together, for example `[GENEC|GENED]`. An intergenic POS appears as
`[]`. All other genes lie to the left or right in genomic coordinates.

The map places each flanking gene at its **nearest body boundary**: its end
when left of POS, its start when right of POS. Overlapping genes share the
POS anchor. Genes sharing another anchor are separated by `|`. This also
handles nested and overlapping gene bodies without pretending to draw their
full lengths.

Dots represent the distance between consecutive anchors and between the
outermost anchors and the window edges. With `--dot-bp 10000`, one dot
represents up to 10 kb: each positive gap is rounded up independently. The
left window edge is clipped at position 1; the right edge is `POS+window`
(annotation contig lengths are not required). Gene labels have variable
width, so this is a schematic, not a fixed-width coordinate plot.

By default a gap is capped at 80 dots; a trailing `~` explicitly marks a
longer gap. Change this with `--max-dots`. Exact distances remain available
in the neighborhood lists. For an empty neighborhood with a 100 bp window,
`--dot-bp 1 --max-dots 3` produces `...~[]...~`.

## INFO fields

`--emit` replaces the default selection; repeated `--emit` options accumulate.
Every group can be requested independently. `--emit all` selects all groups.
The default prefix is `GC_`; `--prefix PC_` is useful when annotating the same
VCF separately with different gene filters.

| `--emit` group | INFO field(s) | Value |
| --- | --- | --- |
| `nearest` | `GC_NEAREST` | Nearest gene-body labels, retaining all ties |
| `distance` | `GC_DISTANCE` | Absolute bp distance to the nearest body; zero inside |
| `nearest-tss` | `GC_TSS` | Genes with the nearest selected TSS, retaining all ties |
| `tss-distance` | `GC_TSS_DISTANCE` | Absolute bp distance to that TSS |
| `neighbors` | `GC_NEIGHBORS` | Genes within `--window`, as `GENEA(-250),GENEB(0),GENEC(+800)` using body distances |
| `neighbors-tss` | `GC_TSS_NEIGHBORS` | Genes with TSSs within `--window`, using the same list syntax |
| `map` | `GC_MAP` | Schematic gene neighborhood with bracketed POS and scaled dots |
| `overlap` | `GC_OVERLAP` | All gene bodies containing POS |
| `ids` | `GC_NEAREST_ID`, `GC_TSS_ID`, `GC_OVERLAP_ID` | Identifiers aligned with the corresponding label lists |
| `biotype` | `GC_BIOTYPE` | Biotypes aligned with the nearest **body** genes |
| `strand` | `GC_STRAND` | Strands aligned with the nearest **body** genes |
| `flanks` | `GC_LEFT`, `GC_LEFT_DISTANCE`, `GC_RIGHT`, `GC_RIGHT_DISTANCE` | Nearest nonoverlapping genes on each genomic side, with absolute bp distances |
| `count` | `GC_N_GENES` | Distinct eligible genes with body distance at most `--window` |
| `overlap-count` | `GC_N_OVERLAP` | Distinct eligible genes whose bodies contain POS |
| `tss-details` | `GC_TSS_DETAIL` | Nearest TSS sites as `GENE\|GENE_ID\|TRANSCRIPT_ID\|TSS\|STRAND\|SOURCE` |
| `promoter` | `GC_PROMOTER`, `GC_PROMOTER_ID`, `GC_PROMOTER_DETAIL` | Promoter-overlapping genes, their IDs, and `GENE\|GENE_ID\|TRANSCRIPT_ID\|TSS\|STRAND\|OFFSET` details |
| `exon` | `GC_EXON` | Genes with an annotated exon containing POS |
| `cds` | `GC_CDS` | Genes with an annotated CDS containing POS |
| `intron` | `GC_INTRON` | Genes with POS between consecutive merged exons of a transcript |
| `feature-details` | `GC_FEATURE` | Overlapping features as `GENE\|GENE_ID\|TRANSCRIPT_ID\|FEATURE\|START\|END` |
| `relations` | `GC_RELATION` | Body-neighborhood labels with gene-oriented `upstream`, `downstream`, `overlap`, or `unknown_strand` relationships |

Distance and count fields use `Number=1,Type=Integer`. The map uses
`Number=1,Type=String`; other fields use `Number=.,Type=String`. The single
nearest distance applies to every tied gene. Missing matches omit the tag.
Known annotation contigs with no neighbors have a count of zero and an empty
variant marker in the map. A contig absent from the annotation leaves all
selected annotations missing, including counts and the map.

Gene symbols are used by default, falling back to identifiers when names are
missing. `--label id` uses identifiers throughout, including neighborhoods
and maps. Distinct identifiers with the same symbol remain separate entries.
Reserved characters in labels are percent-encoded, such as `%2C` for a comma
and `%3B` for a semicolon. Dots in map labels are encoded as `%2E` so they
cannot be confused with distance marks. Identifier versions are preserved.

### Existing annotations and provenance

If a selected INFO tag already exists in the header, the plugin stops unless
`--overwrite` is supplied. Overwrite requires compatible Number/Type and
clears old values when the new query has no match. Unselected tags remain
unchanged. Use a different prefix to retain two annotation sets in the same
record.

A `##genecontext` header records the plugin version, annotation path,
indexed/prepared mode, POS coordinate convention, TSS source and canonical
tag, structural transcript selection, promoter widths, distance/window/map
settings, biotype filters, emitted groups, label/prefix choices, feature
types, and contig-map path. bcftools also records its
command line by default. Its `--no-version` suppresses the host command-line
record, but does not suppress the plugin's scientific settings.

## Performance and limitations

Genecontext implements the bcftools `init()`, `process()`, and `destroy()`
streaming interface. Annotation preparation happens once at startup; each
variant is annotated as it passes through.

For BGZF GTF/GFF3 with an existing index, the plugin attempts to load CSI and
then TBI. It checks that the index uses 1-based chromosome/start/end columns
1/4/5, `#` comments, and no skipped lines. Indexes older than their annotation
are rejected. Invalid or incompatible indexes produce a diagnostic; if no
usable index remains, preparation is automatic. These checks cannot prove
that an index built against a different file with preserved timestamps is
correct; rebuild the index when replacing an annotation.

Without a usable index, the plugin:

1. Streams the plain/gzip/BGZF annotation and retains eligible gene records
   and the transcript/exon/CDS feature types needed by the requested outputs.
2. Sorts them by coordinates with a bounded sort buffer (`--sort-mem 64`, MiB).
3. Merges sorted runs with at most 32 input files open at each merge step.
4. Writes a reduced BGZF annotation table and builds a CSI index using
   HTSlib directly, without invoking external programs.
5. Uses the temporary index throughout the VCF stream and cleans up its private
   scratch directory when finished or on a caught plugin error.

Scratch files are created under `--temp-dir`, otherwise `TMPDIR` or the
system temporary directory. Forced termination, such as SIGKILL, can leave
a `genecontext.*` directory to remove manually. Source annotations and their
indexes are never modified.

The internal table retains coordinates, identifiers, names, biotypes,
feature kinds, parent links, and transcript tags. Coordinates remain 1-based
inclusive. The direct and prepared paths apply the same parser, relationship
resolution, gene filters, and annotation calculations. Child features of
filtered-out genes cannot contribute results. A prepared input is scanned
in full; the direct path validates only queried feature rows.

For transcript-based annotations, the reader retrieves complete annotated
spans of genes intersecting the initial query and passing the gene filters.
This makes transcript validation independent of which parts intersect the
query and supplies flanking exons for intron calculations. Additional genes
outside the original query do not become eligible. Nearest, neighborhood,
and promoter annotations still apply their configured distance limits.

`--sort-mem` is a sort-buffer target, not a total process-memory limit. A
single record, merge buffers, the tabix index, and features in the current
query/model spans require additional memory. Large genes with many isoforms
can require sizable model buffers. Parent relationships are resolved for
each query; a whole-annotation gene/transcript graph is not retained in
memory. Repeated variants at the same contig and POS reuse their
last annotation result. Whole-genome throughput has not been benchmarked.

As with ordinary streaming bcftools commands, a late error can leave a partial
output file. A processing error exits nonzero; it is never interpreted as a
request to silently discard a variant.

These are positional annotations. Functional consequences and causal-gene
predictions require additional interpretation.

## License

[MIT](LICENSE). See [THIRD_PARTY.md](THIRD_PARTY.md) for dependency notices.
