"""Regression coverage for GTF attributes, filtering, and exon aggregation."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import nexons


class ReadGtfTests(unittest.TestCase):
    def parse(self, rows, limit):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.gtf"
            path.write_text("# comment\n" + "".join(
                f'chr1\tx\texon\t{start}\t{end}\t.\t-\t.\t{attrs}\n'
                for start, end, attrs in rows))
            return nexons.read_gtf(path, limit)

    def test_aggregation_and_stable_sort(self):
        attrs = 'gene_id "g"; transcript_id "t"; gene_name "G"; transcript_name "T"; transcript_support_level "2";'
        genes = self.parse([(30, 40, attrs), (10, 20, attrs), (10, 15, attrs)], 2)
        self.assertEqual(genes, {"g": {
            "id": "g", "name": "G", "chrom": "chr1", "strand": "-",
            "start": 10, "end": 40, "transcripts": {"t": {
                "id": "t", "name": "T", "chrom": "chr1", "strand": "-",
                "start": 10, "end": 40, "exons": [[10, 20], [10, 15], [30, 40]]}}}})

    def test_filtering_and_selection_overrides(self):
        base = 'gene_id "g"; transcript_id "t";'
        for tsl in ('', ' transcript_support_level "NA";',
                    ' transcript_support_level "3 (assigned)";'):
            self.assertEqual(self.parse([(10, 20, base + tsl)], 2), {})
            self.assertIn("g", self.parse([(10, 20, base + tsl)], None))
            for tag in ("MANE_Select", "Ensembl_Canonical", "gencode_primary", "gencode_basic"):
                self.assertIn("g", self.parse([(10, 20, base + tsl + f' tag "{tag}";')], 1))
        attrs = base + ' transcript_support_level "2"; transcript_support_level "NA";'
        self.assertIn("g", self.parse([(10, 20, attrs)], 2))

    def test_fallbacks_and_warnings(self):
        genes = self.parse([(10, 20, 'gene_name "g"; transcript_name "t";')], None)
        self.assertEqual(genes["g"]["id"], "g")
        self.assertEqual(genes["g"]["transcripts"]["t"]["id"], "t")
        with patch.object(nexons, "warn") as warning:
            self.assertEqual(self.parse([(10, 20, 'gene_id "g";')], None), {})
            warning.assert_called_once_with("No transcript name or id found for exon at chr1:10-20")
        with patch.object(nexons, "warn") as warning:
            attrs = 'gene_id "g"; transcript_id "t"; transcript_support_level "bad";'
            self.assertEqual(self.parse([(10, 20, attrs)], 2), {})
            warning.assert_called_once_with("Ignoring non-numeric TSL value bad")


if __name__ == "__main__":
    unittest.main()
