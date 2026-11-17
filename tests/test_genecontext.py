"""End-to-end checks against bcftools, using only synthetic annotation data."""
import gzip
import os
from pathlib import Path
import random
import re
import unittest

from helpers import BCFTOOLS, PluginFixture, command, gtf, records


class PluginTests(PluginFixture):
    def test_nearest_bodies_ties_boundaries_and_missing(self):
        out = self.infos(self.run_plugin("--max-distance", "500"))
        for pos in [100, 150, 200]:
            self.assertEqual(out[pos]["GC_NEAREST"], "A")
            self.assertEqual(out[pos]["GC_DISTANCE"], "0")
        self.assertEqual(out[250]["GC_NEAREST"], "A,B")
        self.assertEqual(out[250]["GC_DISTANCE"], "50")
        self.assertEqual(out[475]["GC_NEAREST"], "B,C")
        self.assertNotIn("GC_NEAREST", out[2000])
        self.assertNotIn("GC_DISTANCE", out[2000])
        self.assertTrue(all(x["KEEP"] == "x" for x in out.values()))
        self.assert_clean()

    def test_gene_boundary_tss_on_both_strands(self):
        out = self.infos(self.run_plugin())
        self.assertEqual(out[150]["GC_TSS"], "A")
        self.assertEqual(out[150]["GC_TSS_DISTANCE"], "50")
        self.assertEqual(out[275]["GC_TSS"], "A,C")
        self.assertEqual(out[275]["GC_TSS_DISTANCE"], "175")
        self.assertEqual(out[475]["GC_TSS"], "B,C")
        self.assertEqual(out[475]["GC_TSS_DISTANCE"], "25")
        self.assertEqual(out[1100]["GC_TSS"], "D")
        self.assertEqual(out[1100]["GC_TSS_DISTANCE"], "100")

    def test_neighborhood_lists_and_map(self):
        out = self.infos(self.run_plugin("--emit", "all", "--window", "400", "--dot-bp", "100"))
        self.assertEqual(out[250]["GC_NEIGHBORS"], "A(-50),B(+50),C(+200)")
        self.assertEqual(out[250]["GC_TSS_NEIGHBORS"], "A(-150),C(+200),B(+250)")
        self.assertEqual(out[250]["GC_MAP"], "..A.[].B..C..")
        self.assertEqual(out[475]["GC_MAP"], "..A...[B|C]....")
        self.assertEqual(out[475]["GC_OVERLAP"], "B,C")
        self.assertEqual(out[475]["GC_N_OVERLAP"], "2")
        self.assertEqual(out[475]["GC_N_GENES"], "3")
        self.assertEqual(out[475]["GC_BIOTYPE"], "lncRNA,protein_coding")
        self.assertEqual(out[475]["GC_STRAND"], "-,+")
        self.assertEqual(out[475]["GC_LEFT"], "A")
        self.assertEqual(out[475]["GC_LEFT_DISTANCE"], "275")
        self.assertEqual(out[475]["GC_RIGHT"], "D")
        self.assertEqual(out[475]["GC_RIGHT_DISTANCE"], "425")

    def test_output_groups_are_independent(self):
        for option, tag in [("nearest", "GC_NEAREST"), ("distance", "GC_DISTANCE"),
                            ("nearest-tss", "GC_TSS"), ("tss-distance", "GC_TSS_DISTANCE"),
                            ("map", "GC_MAP"), ("overlap-count", "GC_N_OVERLAP")]:
            with self.subTest(option=option):
                result = self.run_plugin("--emit", option)
                tags = re.findall(r"##INFO=<ID=(GC_[^,]+)", result.stdout.decode())
                self.assertEqual(tags, [tag])

    def test_plain_gzip_and_existing_csi_tbi_agree(self):
        baseline = records(self.run_plugin("--emit", "all").stdout)
        zipped = self.d / "normal.gtf.gz"
        with gzip.open(zipped, "wb") as out:
            out.write(self.genes.read_bytes())
        self.assertEqual(records(self.run_plugin("--emit", "all", annotation=zipped).stdout), baseline)
        for kind in ["csi", "tbi"]:
            indexed = self.indexed_annotation(kind=kind)
            # Direct mode must not need a writable/valid scratch directory.
            result = self.run_plugin("--emit", "all", "--temp-dir", self.d / "absent", annotation=indexed)
            self.assertIn(b"using existing index", result.stderr)
            self.assertEqual(records(result.stdout), baseline)
        self.assert_clean()

    def test_direct_index_ignores_nongene_features(self):
        with self.genes.open("a") as out:
            out.write(gtf(250, 250, "EXON", feature="exon"))
            out.write(gtf(275, 275, "TRANSCRIPT", feature="transcript").rstrip() + ' transcript_id "unused";\n')
        indexed = self.indexed_annotation()
        prepared = self.run_plugin("--emit", "all")
        direct = self.run_plugin("--emit", "all", annotation=indexed)
        self.assertEqual(records(prepared.stdout), records(direct.stdout))
        self.assertEqual(self.infos(direct)[250]["GC_NEAREST"], "A,B")

    def test_gff3_attributes_and_embedded_fasta(self):
        gff = self.d / "genes.gff3"
        gff.write_text("##gff-version 3\n" + "".join(
            f"chr1\ttest\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={name};Name={name};gene_biotype={biotype}\n"
            for start, end, name, strand, biotype in [(100, 200, "A", "+", "protein_coding"),
                                                     (300, 500, "B", "-", "lncRNA"),
                                                     (450, 600, "C", "+", "protein_coding"),
                                                     (900, 1000, "D", "-", "protein_coding")]))
        baseline = records(self.run_plugin("--emit", "all").stdout)
        indexed = self.indexed_annotation(gff)
        self.assertEqual(records(self.run_plugin("--emit", "all", annotation=indexed).stdout), baseline)
        with gff.open("a") as out:
            out.write("##FASTA\n>chr1\nACGTACGT\n")
        self.assertEqual(records(self.run_plugin("--emit", "all", annotation=gff).stdout), baseline)

    def test_biotype_filters_apply_before_nearest_selection(self):
        out = self.infos(self.run_plugin("--protein-coding", "--emit", "all"))
        self.assertEqual(out[300]["GC_NEAREST"], "A")
        self.assertEqual(out[475]["GC_NEAREST"], "C")
        self.assertEqual(out[475]["GC_OVERLAP_ID"], "C")
        indexed = self.indexed_annotation()
        direct = self.infos(self.run_plugin("--protein-coding", "--emit", "all", annotation=indexed))
        self.assertEqual(out, direct)
        out = self.infos(self.run_plugin("--biotypes", "lncRNA"))
        self.assertEqual(out[150]["GC_NEAREST"], "B")
        out = self.infos(self.run_plugin("--exclude-biotypes", "protein_coding"))
        self.assertEqual(out[150]["GC_NEAREST"], "B")

    def test_all_filtered_genes_are_missing_not_an_error(self):
        result = self.run_plugin("--biotypes", "does_not_exist", "--emit", "all")
        self.assertIn(b"no genes pass", result.stderr)
        for _, info in records(result.stdout):
            self.assertNotIn("GC_NEAREST", info)
            self.assertEqual(info["GC_N_GENES"], "0")
        self.assert_clean()

    def test_unknown_strand_and_missing_names(self):
        self.genes.write_text('chr1\tx\tgene\t100\t200\t.\t.\t.\tgene_id "ID.1";\n')
        out = self.infos(self.run_plugin("--emit", "all"))
        self.assertEqual(out[150]["GC_NEAREST"], "ID.1")
        self.assertNotIn("GC_TSS", out[150])
        self.assertNotIn("GC_TSS_DISTANCE", out[150])
        self.assertEqual(out[150]["GC_STRAND"], ".")
        self.assertIn("[ID%2E1]", out[150]["GC_MAP"])

    def test_escape_gtf_gff3_values_and_preserve_distinct_ids(self):
        gff = self.d / "special.gff3"
        gff.write_text("chr1\tx\tgene\t100\t200\t.\t+\t.\tID=id%3B1;Name=A%2CB%28X%29%5BY%5D%7CZ%20Q%25;biotype=protein_coding\n")
        out = self.infos(self.run_plugin("--emit", "all", annotation=gff))
        label = "A%2CB%28X%29%5BY%5D%7CZ%20Q%25"
        self.assertEqual(out[150]["GC_NEAREST"], label)
        self.assertEqual(out[150]["GC_NEAREST_ID"], "id%3B1")
        self.assertIn("[" + label + "]", out[150]["GC_MAP"])
        self.genes.write_text(gtf(100, 200, "SAME", gene_id="one") * 2 + gtf(100, 200, "SAME", gene_id="two"))
        out = self.infos(self.run_plugin("--emit", "all"))
        self.assertEqual(out[150]["GC_NEAREST"], "SAME,SAME")
        self.assertEqual(out[150]["GC_NEAREST_ID"], "one,two")
        out = self.infos(self.run_plugin("--emit", "all", "--label", "id"))
        self.assertEqual(out[150]["GC_NEAREST"], "one,two")

    def test_duplicate_id_parts_use_closest_part_without_filling_gap(self):
        self.genes.write_text(gtf(100, 200, "A", gene_id="one") + gtf(300, 400, "A", gene_id="one"))
        out = self.infos(self.run_plugin("--emit", "all"))
        self.assertEqual(out[250]["GC_DISTANCE"], "50")
        self.assertEqual(out[250]["GC_NEAREST"], "A")
        self.assertNotIn("GC_OVERLAP", out[250])
        self.assertEqual(out[250]["GC_N_GENES"], "1")

    def test_flanks_exclude_every_part_of_an_overlapping_gene(self):
        self.genes.write_text(gtf(1, 50, "LEFT") + gtf(100, 200, "PART") +
                              gtf(300, 400, "PART") + gtf(500, 600, "RIGHT"))
        source = self.write_vcf([100, 150, 200, 250, 300, 350, 400])
        for annotation in [self.genes, self.indexed_annotation(kind="csi"), self.indexed_annotation(kind="tbi")]:
            out = self.infos(self.run_plugin("--emit", "overlap,flanks", source=source, annotation=annotation))
            for pos in [100, 150, 200, 300, 350, 400]:
                with self.subTest(annotation=annotation.name, pos=pos):
                    self.assertEqual(out[pos]["GC_OVERLAP"], "PART")
                    self.assertEqual(out[pos]["GC_LEFT"], "LEFT")
                    self.assertEqual(out[pos]["GC_LEFT_DISTANCE"], str(pos - 50))
                    self.assertEqual(out[pos]["GC_RIGHT"], "RIGHT")
                    self.assertEqual(out[pos]["GC_RIGHT_DISTANCE"], str(500 - pos))
            # An unfilled gap overlaps neither part: both sides can use PART.
            self.assertNotIn("GC_OVERLAP", out[250])
            self.assertEqual(out[250]["GC_LEFT"], "PART")
            self.assertEqual(out[250]["GC_RIGHT"], "PART")
            self.assertEqual(out[250]["GC_LEFT_DISTANCE"], "50")
            self.assertEqual(out[250]["GC_RIGHT_DISTANCE"], "50")
            self.assert_clean()

    def test_multipart_tss_sites_are_unique_despite_interleaved_genes(self):
        for strand, intervals, pos in [("+", [(100, 200), (100, 250), (100, 300)], 100),
                                       ("-", [(100, 400), (150, 400), (200, 400)], 400)]:
            self.genes = self.d / f"multipart-tss-{strand}.gtf"
            self.genes.write_text("".join(gtf(start, end, name, strand, gene_id=gid)
                                         for (start, end), name, gid in zip(intervals, ["A", "B", "A"], ["G", "H", "G"])))
            source = self.write_vcf([pos])
            for annotation in [self.genes, self.indexed_annotation(kind="csi"), self.indexed_annotation(kind="tbi")]:
                with self.subTest(strand=strand, annotation=annotation.name):
                    out = self.infos(self.run_plugin("--emit", "nearest-tss,tss-details,promoter",
                                                    "--promoter-upstream", "0", "--promoter-downstream", "0",
                                                    annotation=annotation, source=source))[pos]
                    self.assertEqual(out["GC_TSS"], "A,B")
                    self.assertEqual(out["GC_TSS_DETAIL"], f"A|G|.|{pos}|{strand}|gene,B|H|.|{pos}|{strand}|gene")
                    self.assertEqual(out["GC_PROMOTER_DETAIL"], f"A|G|.|{pos}|{strand}|0,B|H|.|{pos}|{strand}|0")
                    self.assert_clean()

    def test_zero_and_exact_limits(self):
        out = self.infos(self.run_plugin("--max-distance", "0", "--window", "0", "--emit", "all"))
        self.assertEqual(out[100]["GC_TSS_DISTANCE"], "0")
        self.assertNotIn("GC_TSS", out[150])
        self.assertNotIn("GC_NEAREST", out[250])
        self.assertEqual(out[250]["GC_MAP"], "[]")
        self.assertEqual(out[475]["GC_MAP"], "[B|C]")
        out = self.infos(self.run_plugin("--max-distance", "50", "--window", "50", "--emit", "all"))
        self.assertEqual(out[250]["GC_NEAREST"], "A,B")
        self.assertEqual(out[250]["GC_NEIGHBORS"], "A(-50),B(+50)")
        out = self.infos(self.run_plugin("--max-distance", "49"))
        self.assertNotIn("GC_NEAREST", out[250])

    def test_map_long_gaps_are_explicitly_capped(self):
        out = self.infos(self.run_plugin("--emit", "map", "--window", "100", "--dot-bp", "1", "--max-dots", "3"))
        self.assertEqual(out[2000]["GC_MAP"], "...~[]...~")

    def test_contig_mapping_missing_policy_and_cleanup(self):
        missing = self.write_vcf([100, 150], chrom="chr2", name="missing.vcf")
        result = self.run_plugin("--emit", "all", source=missing)
        self.assertEqual(result.stderr.count(b"annotation has no contig"), 1)
        self.assertEqual([info for _, info in records(result.stdout)], [{"KEEP": "x"}] * 2)
        mapping = self.d / "contigs.tsv"
        mapping.write_text("chr2\tchr1\n")
        result = self.run_plugin("--chr-map", mapping, source=missing)
        self.assertEqual(self.infos(result)[100]["GC_NEAREST"], "A")
        result = self.run_plugin("--missing-contig", "error", source=missing, success=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(b"no contig", result.stderr)
        self.assert_clean()

    def test_incompatible_corrupt_and_stale_indexes_fall_back(self):
        # A zero-based index with the right columns is still incompatible.
        bad = self.indexed_annotation(config=["-0", "-s", "1", "-b", "4", "-e", "5"])
        result = self.run_plugin(annotation=bad)
        self.assertIn(b"incompatible index", result.stderr)
        self.assertIn(b"prepared 4 genes", result.stderr)
        Path(str(bad) + ".csi").write_bytes(b"invalid index")
        result = self.run_plugin(annotation=bad)
        self.assertIn(b"cannot load", result.stderr)
        self.assertEqual(self.infos(result)[250]["GC_NEAREST"], "A,B")
        good = self.indexed_annotation(kind="tbi")
        stat = good.stat()
        os.utime(str(good) + ".tbi", (stat.st_atime - 10, stat.st_mtime - 10))
        result = self.run_plugin(annotation=good)
        self.assertIn(b"older than", result.stderr)
        self.assertEqual(self.infos(result)[250]["GC_NEAREST"], "A,B")
        self.assert_clean()

    def test_streaming_vcf_and_bcf_preserve_format(self):
        source = self.write_vcf([150, 475], samples=True)
        baseline = records(source.read_bytes())
        for typ in ["v", "u"]:
            data = command([BCFTOOLS, "view", "-O" + typ, source]).stdout
            result = self.run_plugin(source="-", data=data)
            annotated = records(result.stdout)
            for (before, _), (after, _) in zip(baseline, annotated):
                self.assertEqual(before[:7], after[:7])
                self.assertEqual(before[8:], after[8:])
            self.assertNotIn(b"##FORMAT=<ID=GC_", result.stdout)
        self.assert_clean()

    def test_output_formats_indexes_and_native_filters(self):
        baseline = records(self.run_plugin().stdout)
        for typ in ["v", "z", "u", "b"]:
            out = self.d / ("out." + typ)
            flags = ["-O" + typ, "-o", out]
            if typ in ["z", "b"]:
                flags += ["-W", "--threads", "2"]
            self.run_plugin(general=flags)
            read = command([BCFTOOLS, "view", "-Ov", out])
            self.assertEqual(records(read.stdout), baseline)
            if typ in ["z", "b"]:
                self.assertTrue(Path(str(out) + ".csi").exists())
                region = command([BCFTOOLS, "view", "-H", "-r", "chr1:250-275", out])
                self.assertEqual([int(c[1]) for c, _ in records(region.stdout)], [250, 275])
        result = self.run_plugin(general=["-i", "POS=250"])
        self.assertEqual(list(self.infos(result)), [250])
        result = self.run_plugin(general=["-t", "chr1:450-475"])
        self.assertEqual(list(self.infos(result)), [450, 475])

    def test_overwrite_compatible_tags_clears_stale_values(self):
        out = self.run_plugin()
        source = self.d / "annotated.vcf"
        source.write_bytes(out.stdout)
        result = self.run_plugin(source=source, success=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(b"already exists", result.stderr)
        result = self.run_plugin("--overwrite", "--max-distance", "0", source=source)
        self.assertNotIn("GC_NEAREST", self.infos(result)[250])
        self.assertEqual(self.infos(result)[150]["GC_NEAREST"], "A")
        result = self.run_plugin("--prefix", "OTHER_", source=source)
        self.assertEqual(self.infos(result)[250]["OTHER_NEAREST"], "A,B")
        self.assertEqual(self.infos(result)[250]["GC_NEAREST"], "A,B")
        source.write_bytes(out.stdout.replace(b"ID=GC_DISTANCE,Number=1,Type=Integer", b"ID=GC_DISTANCE,Number=1,Type=String"))
        result = self.run_plugin("--overwrite", source=source, success=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(b"Incompatible Number/Type", result.stderr)
        self.assert_clean()

    def test_unsorted_input_and_multiple_alleles_use_pos(self):
        source = self.write_vcf([475, 250, 475, 150])
        source.write_text(source.read_text().replace("\tA\tG\t", "\tAAAA\tA,C\t"))
        result = records(self.run_plugin(source=source).stdout)
        self.assertEqual([c[1] for c, _ in result], ["475", "250", "475", "150"])
        self.assertEqual([x["GC_NEAREST"] for _, x in result], ["B,C", "A,B", "B,C", "A"])

    def test_invalid_arguments_and_annotations_fail_cleanly(self):
        for options in [("--window", "-1"), ("--dot-bp", "0"), ("--max-distance", "1.5"),
                        ("--emit", "unknown"), ("--biotypes", "x,,y"), ("--prefix", "1BAD"),
                        ("--protein-coding", "--biotypes", "lncRNA")]:
            with self.subTest(options=options):
                result = self.run_plugin(*options, success=False)
                self.assertNotEqual(result.returncode, 0)
                self.assert_clean()
        for data in [gtf(200, 100, "BAD"), gtf(0, 100, "BAD"),
                     'chr1\tx\tgene\t1\t2\t.\t+\t.\tgene_id "unterminated;\n',
                     'chr1\tx\tgene\t1\t2\t.\t+\t.\tID=broken%ZZ\n',
                     'chr1\tx\tgene\t1\t2\t.\t+\t.\tName=NoID\n',
                     "not a feature\n", gtf(1, 100, "TRANSCRIPT", feature="transcript")]:
            with self.subTest(data=data):
                self.genes.write_text(data)
                result = self.run_plugin(success=False)
                self.assertNotEqual(result.returncode, 0)
                self.assert_clean()

    def test_bad_indexed_gene_causes_nonzero_exit_not_silent_record_loss(self):
        self.genes.write_text('chr1\tx\tgene\t100\t200\t.\t+\t.\tgene_id "broken;\n')
        indexed = self.indexed_annotation()
        result = self.run_plugin(annotation=indexed, success=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(b"Unterminated GTF", result.stderr)
        self.assert_clean()

    def test_csi_long_coordinates(self):
        self.genes.write_text(gtf(700000000, 700000100, "LONG", "-"))
        source = self.write_vcf([700000150])
        for annotation in [self.genes, self.indexed_annotation()]:
            result = self.infos(self.run_plugin(source=source, annotation=annotation))[700000150]
            self.assertEqual(result["GC_DISTANCE"], "50")
            self.assertEqual(result["GC_TSS_DISTANCE"], "50")

    def test_multilevel_sort_merge_and_scratch_failure(self):
        # Force more than 32 runs to exercise bounded fan-in, with a small
        # reference whose output is independently predictable.
        with self.genes.open("w") as out:
            for i in reversed(range(1, 7001)):
                out.write(gtf(i * 10, i * 10 + 5, "X" * 5000, gene_id=f"g{i}"))
        source = self.write_vcf([10, 16, 70005])
        result = self.run_plugin("--sort-mem", "1", "--emit", "distance", source=source)
        self.assertGreater(int(re.search(rb"\((\d+) sort runs\)", result.stderr)[1]), 32)
        self.assertEqual([info["GC_DISTANCE"] for _, info in records(result.stdout)], ["0", "1", "0"])
        self.assert_clean()
        result = self.run_plugin("--temp-dir", self.d / "nonexistent", success=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(b"Cannot create temporary directory", result.stderr)

    def test_corrupt_gzip_annotation_fails_and_cleans_up(self):
        path = self.d / "broken.gtf.gz"
        data = gzip.compress(self.genes.read_bytes())
        path.write_bytes(data[:len(data) // 2])
        result = self.run_plugin(annotation=path, success=False)
        self.assertNotEqual(result.returncode, 0)
        self.assert_clean()

    def test_spilled_sort_and_randomized_distance_oracle(self):
        rng = random.Random(42)
        genes = []
        for i in range(8000):
            start = rng.randrange(1, 1000000)
            genes.append((start, start + rng.randrange(1, 1500), f"g{i:05}", rng.choice("+-")))
        self.genes.write_text("".join(gtf(start, end, name, strand) for start, end, name, strand in genes))
        positions = [rng.randrange(1, 1000000) for _ in range(40)]
        source = self.write_vcf(positions)
        small = self.run_plugin("--sort-mem", "1", "--max-distance", "100", source=source)
        large = self.run_plugin("--sort-mem", "64", "--max-distance", "100", source=source)
        self.assertEqual(records(small.stdout), records(large.stdout))
        self.assertRegex(small.stderr.decode(), r"\([2-9][0-9]* sort runs\)")
        for columns, info in records(small.stdout):
            pos = int(columns[1])
            for metric, name_tag, dist_tag in [("body", "GC_NEAREST", "GC_DISTANCE"),
                                              ("tss", "GC_TSS", "GC_TSS_DISTANCE")]:
                distances = [(max(start-pos, pos-end, 0) if metric == "body" else abs(pos-(start if strand == "+" else end)),
                              start, end, name) for start, end, name, strand in genes]
                distance = min(x[0] for x in distances)
                if distance <= 100:
                    expected = ",".join(x[3] for x in sorted(distances) if x[0] == distance)
                    self.assertEqual(info[name_tag], expected)
                    self.assertEqual(info[dist_tag], str(distance))
                else:
                    self.assertNotIn(name_tag, info)
        self.assert_clean()


if __name__ == "__main__":
    unittest.main()
