"""Shared synthetic fixtures and bcftools helpers for integration tests."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


BCFTOOLS = os.environ.get("BCFTOOLS", "bcftools")
BGZIP = os.environ.get("BGZIP", "bgzip")
TABIX = os.environ.get("TABIX", "tabix")
PROJECT_ROOT = Path(__file__).resolve().parents[1]
PLUGIN = Path(os.environ.get("GC_PLUGIN", PROJECT_ROOT / "build/genecontext.so")).resolve()


def command(args, data=None, success=True):
    result = subprocess.run([str(x) for x in args], input=data, capture_output=True, timeout=120)
    if success and result.returncode:
        raise AssertionError(f"Command failed: {args}\n{result.stderr.decode(errors='replace')}")
    return result


def records(data):
    if isinstance(data, bytes):
        data = data.decode()
    result = []
    for line in data.splitlines():
        if line.startswith("#"):
            continue
        columns = line.split("\t")
        info = dict(x.split("=", 1) if "=" in x else (x, True)
                    for x in columns[7].split(";") if x != ".")
        result.append((columns, info))
    return result


def gtf(start, end, name, strand="+", biotype="protein_coding", chrom="chr1", gene_id=None, feature="gene"):
    return (f'{chrom}\ttest\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t'
            f'gene_id "{gene_id or name}"; gene_name "{name}"; gene_type "{biotype}";\n')


class PluginFixture(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="genecontext-tests-")
        self.addCleanup(self.temp.cleanup)
        self.d = Path(self.temp.name)
        self.scratch = self.d / "scratch"
        self.scratch.mkdir()
        self.genes = self.d / "genes.gtf"
        self.genes.write_text(gtf(900, 1000, "D", "-") + gtf(300, 500, "B", "-", "lncRNA") +
                              gtf(100, 200, "A") + gtf(450, 600, "C"))
        self.vcf = self.write_vcf([50, 100, 150, 200, 250, 275, 300, 450, 475, 600, 700, 1100, 2000])

    def write_vcf(self, positions, name="input.vcf", chrom="chr1", extra="", samples=False):
        header = ('##fileformat=VCFv4.3\n##contig=<ID=chr1,length=1000000000>\n'
                  '##contig=<ID=chr2,length=1000000000>\n'
                  '##INFO=<ID=KEEP,Number=1,Type=String,Description="Preserve me">\n' + extra)
        if samples:
            header += ('##FORMAT=<ID=LP,Number=1,Type=Float,Description="LP">\n'
                       '##FORMAT=<ID=LF_START,Number=1,Type=Integer,Description="Locus start">\n')
        header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        header += "\tFORMAT\tTRAIT1\tTRAIT2\n" if samples else "\n"
        for i, pos in enumerate(positions):
            header += f"{chrom}\t{pos}\tv{i}\tA\tG\t.\tPASS\tKEEP=x"
            header += "\tLP:LF_START\t9:50\t.:.\n" if samples else "\n"
        path = self.d / name
        path.write_text(header)
        return path

    def run_plugin(self, *options, annotation=None, source=None, general=(), data=None, success=True):
        return command([BCFTOOLS, "+" + str(PLUGIN), *general, source or self.vcf, "--",
                        "-a", annotation or self.genes, "--temp-dir", self.scratch, *options],
                       data=data, success=success)

    def infos(self, result):
        return {int(cols[1]): info for cols, info in records(result.stdout)}

    def indexed_annotation(self, source=None, kind="csi", config=None):
        source = source or self.genes
        lines = source.read_text().splitlines()
        header = [x for x in lines if x.startswith("#") and not x.startswith("##FASTA")]
        body = sorted((x for x in lines if not x.startswith("#")),
                      key=lambda x: (x.split("\t")[0], int(x.split("\t")[3]), int(x.split("\t")[4]), x))
        path = self.d / (source.name + "." + kind + ".gz")
        data = ("\n".join(header + body) + "\n").encode()
        path.write_bytes(command([BGZIP, "-c"], data=data).stdout)
        flags = ["-C"] if kind == "csi" else []
        command([TABIX, *flags, *(config or ["-p", "gff"]), path])
        return path

    def assert_clean(self):
        self.assertEqual(list(self.scratch.iterdir()), [])
