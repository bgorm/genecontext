"""Transcript-parent resolution and independent coordinate expectations."""
import random

from helpers import PluginFixture, records


GENES = [
    (1000, 3000, "G1", "A", "+", "protein_coding"),
    (4000, 6000, "G2", "B", "-", "lncRNA"),
    (7000, 8000, "G3", "C", "+", "protein_coding"),
    (9000, 10000, "G4", "U", ".", "protein_coding"),
]
TRANSCRIPTS = [
    (1200, 2800, "T1", "G1", "+", ["basic", "Ensembl_canonical"]),
    (1000, 3000, "T2", "G1", "+", ["basic", "MANE_Select"]),
    (4200, 5700, "T3", "G2", "-", ["basic", "Ensembl_canonical"]),
    (7100, 7900, "T4", "G3", "+", ["basic"]),
    (9200, 9800, "T5", "G4", ".", ["basic", "Ensembl_canonical"]),
]
SEGMENTS = [
    (1200, 1400, "T1", "G1", "exon"), (2000, 2200, "T1", "G1", "exon"),
    (2600, 2800, "T1", "G1", "exon"), (1300, 1400, "T1", "G1", "CDS"),
    (2000, 2100, "T1", "G1", "CDS"), (2650, 2700, "T1", "G1", "CDS"),
    (1000, 1100, "T2", "G1", "exon"), (1350, 2050, "T2", "G1", "exon"),
    (2900, 3000, "T2", "G1", "exon"), (1400, 2000, "T2", "G1", "CDS"),
    (4200, 4500, "T3", "G2", "exon"), (5000, 5200, "T3", "G2", "exon"),
    (5500, 5700, "T3", "G2", "exon"), (7100, 7200, "T4", "G3", "exon"),
    (7800, 7900, "T4", "G3", "exon"), (9200, 9400, "T5", "G4", "exon"),
    (9600, 9800, "T5", "G4", "exon"),
]


class TranscriptTests(PluginFixture):
    def setUp(self):
        super().setUp()
        self.write_models()
        self.vcf = self.write_vcf([999, 1000, 1100, 1150, 1200, 1250, 1251, 1300, 1400,
                                   1401, 1500, 1999, 2000, 2201, 2801, 2850, 3000, 3450,
                                   4100, 4501, 5649, 5650, 5700, 5900, 5901, 6100,
                                   7050, 7500, 7950, 8900, 9300, 9500])

    def write_models(self, fmt="gtf", genes=GENES, transcripts=TRANSCRIPTS, segments=SEGMENTS):
        lines = []
        strands = {x[2]: x[4] for x in genes}
        for start, end, gid, name, strand, biotype in genes:
            attrs = (f'gene_id "{gid}"; gene_name "{name}"; gene_type "{biotype}";' if fmt == "gtf" else
                     f"ID=gene:{gid};gene_id={gid};Name={name};biotype={biotype}")
            lines.append(f"chr1\ttest\tgene\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n")
        for start, end, tid, gid, strand, tags in transcripts:
            if fmt == "gtf":
                attrs = f'gene_id "{gid}"; transcript_id "{tid}"; ' + " ".join(f'tag "{tag}";' for tag in tags)
                kind = "transcript"
            else:
                attrs = f"ID=transcript:{tid};transcript_id={tid};Parent=gene:{gid};tag=" + ",".join(tags)
                kind = "mRNA" if gid == "G1" else "ncRNA"
            lines.append(f"chr1\ttest\t{kind}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n")
        for start, end, tid, gid, kind in segments:
            attrs = f'gene_id "{gid}"; transcript_id "{tid}";' if fmt == "gtf" else f"Parent=transcript:{tid}"
            lines.append(f"chr1\ttest\t{kind}\t{start}\t{end}\t.\t{strands[gid]}\t.\t{attrs}\n")
        random.Random(12).shuffle(lines)
        path = self.d / ("model." + fmt)
        path.write_text("".join(lines))
        if fmt == "gtf":
            self.genes = path
        return path

    def test_canonical_and_gene_tss_differ_on_both_strands(self):
        gene = self.infos(self.run_plugin("--emit", "all"))
        canonical = self.infos(self.run_plugin("--tss-source", "canonical", "--emit", "all"))
        self.assertEqual(gene[1250]["GC_TSS_DISTANCE"], "250")
        self.assertEqual(canonical[1250]["GC_TSS_DISTANCE"], "50")
        self.assertEqual(canonical[1250]["GC_TSS_DETAIL"], "A|G1|T1|1200|+|canonical")
        self.assertEqual(gene[4100]["GC_TSS_DISTANCE"], "1900")
        self.assertEqual(canonical[4100]["GC_TSS_DISTANCE"], "1600")
        self.assertEqual(canonical[4100]["GC_TSS_DETAIL"], "B|G2|T3|5700|-|canonical")
        self.assertEqual(canonical[3450]["GC_TSS"], "A,B")
        self.assertEqual(canonical[3450]["GC_TSS_DISTANCE"], "2250")
        self.assertEqual(canonical[3450]["GC_TSS_DETAIL"], "A|G1|T1|1200|+|canonical,B|G2|T3|5700|-|canonical")
        self.assertEqual(gene[1250]["GC_NEAREST"], canonical[1250]["GC_NEAREST"])

    def test_canonical_tag_override_and_no_implicit_fallback(self):
        out = self.infos(self.run_plugin("--tss-source", "canonical", "--canonical-tag", "MANE_Select",
                                        "--emit", "nearest-tss,tss-distance,tss-details", "--max-distance", "1000"))
        self.assertEqual(out[1250]["GC_TSS_DETAIL"], "A|G1|T2|1000|+|canonical")
        self.assertNotIn("GC_TSS", out[5700])
        out = self.infos(self.run_plugin("--tss-source", "canonical", "--max-distance", "1000", "--emit", "all"))
        self.assertEqual(out[7500]["GC_NEAREST"], "C")
        self.assertNotIn("GC_TSS", out[7500])
        self.assertNotIn("GC_TSS", out[9500])
        self.assertEqual(out[9500]["GC_INTRON"], "U")
        self.assertNotIn("GC_PROMOTER", out[9500])

    def test_promoter_windows_are_strand_aware_inclusive_and_independent(self):
        out = self.infos(self.run_plugin("--tss-source", "canonical", "--emit", "promoter",
                                        "--promoter-upstream", "200", "--promoter-downstream", "50",
                                        "--max-distance", "0", "--window", "0"))
        for pos in [1000, 1100, 1150, 1200, 1250]:
            self.assertEqual(out[pos]["GC_PROMOTER"], "A")
        for pos in [5650, 5700, 5900]:
            self.assertEqual(out[pos]["GC_PROMOTER"], "B")
        for pos in [999, 1251, 5649, 5901]:
            self.assertNotIn("GC_PROMOTER", out[pos])
        self.assertEqual(out[1000]["GC_PROMOTER_DETAIL"], "A|G1|T1|1200|+|-200")
        self.assertEqual(out[5650]["GC_PROMOTER_DETAIL"], "B|G2|T3|5700|-|50")
        self.assertEqual(out[5900]["GC_PROMOTER_DETAIL"], "B|G2|T3|5700|-|-200")
        self.assertEqual(out[5900]["GC_PROMOTER_ID"], "G2")

    def test_gene_boundary_and_zero_width_promoters(self):
        out = self.infos(self.run_plugin("--emit", "promoter", "--promoter-upstream", "0", "--promoter-downstream", "0"))
        self.assertEqual(out[1000]["GC_PROMOTER_DETAIL"], "A|G1|.|1000|+|0")
        self.assertNotIn("GC_PROMOTER", out[1200])
        out = self.infos(self.run_plugin("--tss-source", "canonical", "--emit", "promoter",
                                        "--promoter-upstream", "0", "--promoter-downstream", "0"))
        self.assertEqual(out[1200]["GC_PROMOTER"], "A")
        self.assertEqual(out[5700]["GC_PROMOTER"], "B")
        self.assertNotIn("GC_PROMOTER", out[1000])

    def test_exon_cds_and_intron_preserve_alternative_transcript_context(self):
        out = self.infos(self.run_plugin("--emit", "exon,cds,intron,feature-details"))
        self.assertEqual(out[1500]["GC_EXON"], "A")
        self.assertEqual(out[1500]["GC_CDS"], "A")
        self.assertEqual(out[1500]["GC_INTRON"], "A")
        self.assertEqual(set(out[1500]["GC_FEATURE"].split(",")), {
            "A|G1|T1|intron|1401|1999", "A|G1|T2|exon|1350|2050", "A|G1|T2|CDS|1400|2000"})
        self.assertIn("A|G1|T1|exon|1200|1400", out[1400]["GC_FEATURE"])
        self.assertIn("A|G1|T1|intron|1401|1999", out[1401]["GC_FEATURE"])
        self.assertIn("A|G1|T1|intron|1401|1999", out[1999]["GC_FEATURE"])
        self.assertIn("A|G1|T1|exon|2000|2200", out[2000]["GC_FEATURE"])
        self.assertEqual(out[4501]["GC_FEATURE"], "B|G2|T3|intron|4501|4999")

    def test_canonical_feature_subset_and_transcript_tails(self):
        out = self.infos(self.run_plugin("--emit", "exon,cds,intron,feature-details", "--feature-transcripts", "canonical"))
        self.assertEqual(out[1500]["GC_FEATURE"], "A|G1|T1|intron|1401|1999")
        self.assertNotIn("GC_EXON", out[1500])
        self.assertNotIn("GC_CDS", out[1500])
        for pos in [1100, 1150, 2801, 2850, 7050, 7500, 7950]:
            self.assertNotIn("GC_FEATURE", out[pos])
        all_tx = self.infos(self.run_plugin("--emit", "intron"))
        self.assertEqual(all_tx[7500]["GC_INTRON"], "C")
        self.assertNotIn("GC_INTRON", all_tx[7050])
        self.assertNotIn("GC_INTRON", all_tx[7950])

    def test_gtf_gff3_and_both_index_types_agree_for_all_extended_annotations(self):
        options = ["--emit", "all", "--tss-source", "canonical", "--window", "1000"]
        baseline = records(self.run_plugin(*options).stdout)
        for source in [self.genes, self.write_models("gff3")]:
            self.assertEqual(records(self.run_plugin(*options, annotation=source).stdout), baseline)
            for kind in ["csi", "tbi"]:
                indexed = self.indexed_annotation(source, kind)
                result = self.run_plugin(*options, "--temp-dir", self.d / "absent", annotation=indexed)
                self.assertIn(b"using existing index", result.stderr)
                self.assertEqual(records(result.stdout), baseline)
        self.assert_clean()

    def test_gencode_unquoted_attributes_in_gene_and_child_records(self):
        options = ["--emit", "all", "--tss-source", "canonical", "--window", "1000"]
        baseline = records(self.run_plugin(*options).stdout)
        lines = []
        for line in self.genes.read_text().splitlines():
            columns = line.split("\t")
            # GENCODE supplies bare numeric level/exon_number values alongside
            # quoted identifiers and repeated transcript tags.
            columns[8] = "level 2; " + columns[8]
            if columns[2] in ("exon", "CDS"):
                columns[8] += " exon_number 1;"
            lines.append("\t".join(columns))
        self.genes.write_text("\n".join(lines) + "\n")
        for annotation in [self.genes, self.indexed_annotation()]:
            with self.subTest(annotation=annotation.name):
                result = self.run_plugin(*options, annotation=annotation)
                self.assertEqual(records(result.stdout), baseline)
                out = self.infos(result)
                self.assertEqual(out[1250]["GC_TSS_DETAIL"], "A|G1|T1|1200|+|canonical")
                self.assertEqual(out[1500]["GC_INTRON"], "A")
        self.assert_clean()

    def test_unquoted_gtf_values_still_require_values_and_separators(self):
        for index, suffix in enumerate(['level;', 'level ;', 'level', 'level 2 gene_name "B";',
                                        'level 2";', 'level 2\\3;']):
            with self.subTest(suffix=suffix):
                source = self.d / f"invalid-{index}.gtf"
                source.write_text(
                    'chr1\ttest\tgene\t1000\t3000\t.\t+\t.\tgene_id "G1"; ' + suffix + '\n')
                for annotation in [source, self.indexed_annotation(source)]:
                    result = self.run_plugin(annotation=annotation, success=False)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn(b"GTF", result.stderr)
                    self.assert_clean()

    def test_intron_context_fetches_exons_outside_zero_length_query(self):
        # --emit intron needs only POS for discovery, but both flanking exons
        # are far away. The indexed path must retrieve the complete model.
        options = ["--emit", "intron,feature-details", "--max-distance", "0", "--window", "0",
                   "--feature-transcripts", "canonical"]
        for annotation in [self.genes, self.indexed_annotation()]:
            out = self.infos(self.run_plugin(*options, annotation=annotation))
            self.assertEqual(out[1500]["GC_FEATURE"], "A|G1|T1|intron|1401|1999")
            self.assertEqual(out[4501]["GC_FEATURE"], "B|G2|T3|intron|4501|4999")

    def test_child_biotypes_do_not_override_gene_filter(self):
        options = ["--emit", "all", "--tss-source", "canonical", "--protein-coding"]
        for annotation in [self.genes, self.indexed_annotation()]:
            out = self.infos(self.run_plugin(*options, annotation=annotation))
            self.assertNotIn("GC_INTRON", out[4501])
            self.assertNotIn("GC_EXON", out[5700])
            self.assertNotIn("GC_PROMOTER", out[5700])
            self.assertEqual(out[4100]["GC_TSS"], "A")

    def test_parent_lists_and_percent_escaped_transcript_identifiers(self):
        source = self.write_models("gff3")
        text = source.read_text().replace("transcript:T1", "transcript:T%2C1").replace("transcript_id=T1;", "transcript_id=T%2C1;")
        text += "chr1\tx\texon\t1500\t1550\t.\t+\t.\tParent=transcript:T%2C1,transcript:T2\n"
        source.write_text(text)
        for annotation in [source, self.indexed_annotation(source)]:
            out = self.infos(self.run_plugin("--emit", "all", "--tss-source", "canonical", annotation=annotation))
            self.assertEqual(out[1250]["GC_TSS_DETAIL"], "A|G1|T%2C1|1200|+|canonical")
            details = out[1500]["GC_FEATURE"].split(",")
            self.assertIn("A|G1|T%2C1|exon|1500|1550", details)
            self.assertIn("A|G1|T2|exon|1500|1550", details)
            self.assertNotIn("GC_INTRON", out[1500])

    def test_missing_transcript_rows_can_use_explicit_gtf_exon_links_for_introns(self):
        self.write_models(transcripts=[])
        out = self.infos(self.run_plugin("--emit", "all", "--tss-source", "canonical"))
        self.assertNotIn("GC_TSS", out[1500])
        self.assertEqual(out[1500]["GC_INTRON"], "A")
        self.assertIn("A|G1|T1|intron|1401|1999", out[1500]["GC_FEATURE"])
        out = self.infos(self.run_plugin("--emit", "feature-details", "--feature-transcripts", "canonical"))
        self.assertNotIn("GC_FEATURE", out[1500])

    def test_gene_oriented_relations(self):
        out = self.infos(self.run_plugin("--emit", "relations", "--window", "1000"))
        self.assertEqual(out[6100]["GC_RELATION"], "B(upstream),C(upstream)")
        self.assertIn("U(unknown_strand)", out[8900]["GC_RELATION"])
        self.assertEqual(out[9300]["GC_RELATION"], "U(overlap)")

    def test_overlapping_and_adjacent_exons_do_not_create_false_introns(self):
        segments = [(1200, 1500, "T1", "G1", "exon"), (1400, 1600, "T1", "G1", "exon"),
                    (1601, 1800, "T1", "G1", "exon"), (1900, 2000, "T1", "G1", "exon")]
        self.write_models(segments=segments)
        source = self.write_vcf([1500, 1600, 1601, 1800, 1801, 1899, 1900, 2100])
        out = self.infos(self.run_plugin("--emit", "exon,intron,feature-details", source=source))
        for pos in [1500, 1600, 1601, 1800, 1900, 2100]:
            self.assertNotIn("GC_INTRON", out[pos])
        for pos in [1801, 1899]:
            self.assertEqual(out[pos]["GC_FEATURE"], "A|G1|T1|intron|1801|1899")

    def test_multiple_tagged_transcripts_keep_all_tied_sites(self):
        transcripts = list(TRANSCRIPTS) + [(1300, 2800, "T6", "G1", "+", ["Ensembl_canonical"])]
        self.write_models(transcripts=transcripts)
        out = self.infos(self.run_plugin("--emit", "all", "--tss-source", "canonical"))
        self.assertEqual(out[1250]["GC_TSS"], "A")
        self.assertEqual(out[1250]["GC_TSS_DETAIL"], "A|G1|T1|1200|+|canonical,A|G1|T6|1300|+|canonical")
        self.assertEqual(out[1250]["GC_TSS_NEIGHBORS"].split(",")[0], "A(-50)")

    def test_new_options_are_validated_and_fail_cleanly(self):
        for options in [("--tss-source", "longest"), ("--canonical-tag", ""),
                        ("--promoter-upstream", "-1"), ("--promoter-downstream", "1.5"),
                        ("--feature-transcripts", "unknown"), ("--transcript-features", "gene", "--emit", "all")]:
            with self.subTest(options=options):
                result = self.run_plugin(*options, success=False)
                self.assertNotEqual(result.returncode, 0)
                self.assert_clean()

    def test_extended_info_survives_bcf_stream_and_overwrite(self):
        first = self.run_plugin("--emit", "all", "--tss-source", "canonical", general=["-Ou"])
        result = self.run_plugin("--emit", "all", "--overwrite", "--tss-source", "gene", source="-", data=first.stdout)
        self.assertEqual(self.infos(result)[1250]["GC_TSS_DETAIL"], "A|G1|.|1000|+|gene")
        self.assertIn(b'TSS="gene_boundary"', result.stdout)
        self.assertIn(b'CanonicalTag="Ensembl_canonical"', result.stdout)
        self.assert_clean()

    def test_conflicting_transcript_models_fail_with_cleanup(self):
        # The same transcript ID must not silently choose whichever coordinates
        # happened to appear first in the source or in a tabix query.
        transcripts = list(TRANSCRIPTS) + [(1300, 2800, "T1", "G1", "+", ["Ensembl_canonical"])]
        self.write_models(transcripts=transcripts)
        for annotation in [self.genes, self.indexed_annotation()]:
            result = self.run_plugin("--emit", "all", "--tss-source", "canonical", annotation=annotation, success=False)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(b"Conflicting or multipart transcript coordinates", result.stderr)
            self.assert_clean()

    def test_child_without_transcript_is_not_used_to_invent_introns(self):
        source = self.d / "direct.gff3"
        source.write_text("chr1\tx\tgene\t1000\t3000\t.\t+\t.\tID=gene:G1;gene_id=G1;Name=A\n"
                          "chr1\tx\texon\t1200\t1400\t.\t+\t.\tParent=gene:G1\n"
                          "chr1\tx\texon\t2000\t2200\t.\t+\t.\tParent=gene:G1\n"
                          "chr1\tx\tCDS\t1300\t1400\t.\t+\t0\tParent=gene:G1\n")
        for annotation in [source, self.indexed_annotation(source)]:
            out = self.infos(self.run_plugin("--emit", "exon,cds,intron,feature-details", annotation=annotation))
            self.assertEqual(out[1300]["GC_EXON"], "A")
            self.assertEqual(out[1300]["GC_CDS"], "A")
            self.assertIn("A|G1|.|exon|1200|1400", out[1300]["GC_FEATURE"])
            self.assertNotIn("GC_INTRON", out[1500])
