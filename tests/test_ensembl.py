"""Regression checks using unmodified Ensembl 116 gene models."""
import csv
import gzip
from pathlib import Path

from helpers import PluginFixture, records


FIXTURES = Path(__file__).resolve().parent / "fixtures"

# Expected intervals were checked against the provider rows, independently of
# plugin output. Introns are inclusive gaps between the specified exons.
# The final flag says whether the transcript is tagged Ensembl_canonical.
STRUCTURES = [
    ("1", 12010, "DDX11L1|ENSG00000223972|ENST00000450305|exon|12010|12057", True),
    ("1", 12100, "DDX11L1|ENSG00000223972|ENST00000450305|intron|12058|12178", True),
    ("1", 14900, "WASH7P|ENSG00000227232|ENST00000488147|intron|14830|14969", True),
    ("1", 15000, "WASH7P|ENSG00000227232|ENST00000488147|exon|14970|15038", True),
    ("10", 22475455, "ENSG00000285254|ENSG00000285254|ENST00000642887|exon|22475455|22478026", False),
    ("11", 5225000, "ENSG00000285498|ENSG00000285498|ENST00000644706|exon|5224309|5225461", True),
    ("2", 88857683, "IGKC|ENSG00000211592|ENST00000390237|CDS|88857361|88857683", True),
    ("22", 22030934, "IGLV4-69|ENSG00000211637|ENST00000390282|exon|22030934|22031042", True),
    ("22", 22030994, "IGLV4-69|ENSG00000211637|ENST00000390282|CDS|22030994|22031042", True),
    ("22", 22031100, "IGLV4-69|ENSG00000211637|ENST00000390282|intron|22031043|22031162", True),
    ("22", 22893692, "IGLJ1|ENSG00000211674|ENST00000390320|CDS|22893692|22893818", True),
    ("22", 22895375, "IGLC1|ENSG00000211675|ENST00000390321|CDS|22895375|22895694", True),
    ("7", 142795705, "ENSG00000292271|ENSG00000292271|ENST00000710614|CDS|142795705|142795720", True),
]


class EnsemblTests(PluginFixture):
    def setUp(self):
        super().setUp()
        self.annotation = FIXTURES / "ensembl116.gff3"
        with (FIXTURES / "ensembl116-canonical.tsv").open() as source:
            self.sites = list(csv.DictReader(source, delimiter="\t"))
        self.compressed = self.d / "ensembl116.gff3.gz"
        self.compressed.write_bytes(gzip.compress(self.annotation.read_bytes(), mtime=0))
        self.inputs = [self.annotation, self.compressed,
                       self.indexed_annotation(self.annotation, "csi"),
                       self.indexed_annotation(self.annotation, "tbi")]

    def write_variants(self, positions):
        header = "##fileformat=VCFv4.3\n"
        for line in self.annotation.read_text().splitlines():
            if line.startswith("##sequence-region"):
                _, chrom, _, length = line.split()
                header += f"##contig=<ID={chrom},length={length}>\n"
        header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
        for chrom, pos in sorted(set(positions)):
            header += f"{chrom}\t{pos}\tprobe_{chrom}_{pos}\tA\tG\t.\tPASS\t.\n"
        path = self.d / "ensembl-probes.vcf"
        path.write_text(header)
        return path

    def annotate(self, annotation, variants, *options):
        result = self.run_plugin(*options, annotation=annotation, source=variants)
        indexed = annotation.name.endswith((".csi.gz", ".tbi.gz"))
        self.assertIn(b"using existing index" if indexed else b"prepared", result.stderr)
        self.assert_clean()
        return records(result.stdout)

    def test_canonical_starts_promoters_and_exons(self):
        variants = self.write_variants((site["chrom"], int(site["pos"])) for site in self.sites)
        for annotation in self.inputs:
            result = self.annotate(annotation, variants, "--tss-source", "canonical",
                                   "--feature-transcripts", "canonical", "--emit",
                                   "nearest-tss,tss-distance,tss-details,promoter,exon",
                                   "--max-distance", "0", "--window", "0",
                                   "--promoter-upstream", "0", "--promoter-downstream", "0")
            infos = {(cols[0], cols[1]): info for cols, info in result}
            for site in self.sites:
                with self.subTest(annotation=annotation.name, gene=site["gene_id"]):
                    info = infos[(site["chrom"], site["pos"])]
                    detail = "|".join(site[key] for key in ("gene", "gene_id", "transcript_id", "pos", "strand"))
                    self.assertIn(detail + "|canonical", info.get("GC_TSS_DETAIL", "").split(","))
                    self.assertEqual(info.get("GC_TSS_DISTANCE"), "0")
                    self.assertIn(detail + "|0", info.get("GC_PROMOTER_DETAIL", "").split(","))
                    self.assertIn(site["gene"], info.get("GC_EXON", "").split(","))

    def test_structural_parent_links_and_canonical_selection(self):
        variants = self.write_variants((chrom, pos) for chrom, pos, _, _ in STRUCTURES)
        for annotation in self.inputs:
            for selection in ("all", "canonical"):
                result = self.annotate(annotation, variants, "--emit", "exon,cds,intron,feature-details",
                                       "--feature-transcripts", selection, "--max-distance", "0", "--window", "0")
                infos = {(cols[0], int(cols[1])): info for cols, info in result}
                for chrom, pos, detail, canonical in STRUCTURES:
                    with self.subTest(annotation=annotation.name, selection=selection, detail=detail):
                        info = infos[(chrom, pos)]
                        details = info.get("GC_FEATURE", "").split(",")
                        if selection == "canonical" and not canonical:
                            self.assertNotIn(detail, details)
                        else:
                            self.assertIn(detail, details)
                            gene, _, _, feature, _, _ = detail.split("|")
                            self.assertIn(gene, info.get("GC_" + feature.upper(), "").split(","))

    def test_all_outputs_match_between_preparation_and_existing_indexes(self):
        positions = {(site["chrom"], int(site["pos"]) + offset)
                     for site in self.sites for offset in (-1, 0, 1)}
        for chrom, pos, detail, _ in STRUCTURES:
            positions.add((chrom, pos))
            for boundary in map(int, detail.split("|")[-2:]):
                positions.update((chrom, boundary + offset) for offset in (-1, 0, 1))
        variants = self.write_variants(positions)
        for tss, transcripts in (("gene", "all"), ("canonical", "canonical")):
            baseline = None
            for annotation in self.inputs:
                with self.subTest(tss=tss, transcripts=transcripts, annotation=annotation.name):
                    result = self.annotate(annotation, variants, "--emit", "all", "--tss-source", tss,
                                           "--feature-transcripts", transcripts, "--max-distance", "1000", "--window", "500")
                    if baseline is None:
                        baseline = result
                    self.assertEqual(result, baseline)
