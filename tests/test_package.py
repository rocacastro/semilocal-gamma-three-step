"""Consistency checks for the standalone supplementary package."""
from __future__ import annotations
import json
from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'supplementary_material_gamma.tex'


class PackageTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = SOURCE.read_text(encoding='utf-8')

    def test_six_supplementary_sections(self):
        self.assertEqual(len(re.findall(r'\\section\{', self.text)), 6)

    def test_labels_resolve(self):
        labels = re.findall(r'\\label\{([^}]+)\}', self.text)
        self.assertEqual(len(labels), len(set(labels)))
        used = re.findall(r'\\(?:eqref|ref)\{([^}]+)\}', self.text)
        self.assertFalse(set(used)-set(labels))

    def test_citations_resolve(self):
        entries = re.findall(r'\\bibitem\{([^}]+)\}', self.text)
        self.assertEqual(len(entries), len(set(entries)))
        keys = []
        for group in re.findall(r'\\cite(?:\[[^\]]*\])?\{([^}]+)\}', self.text):
            keys.extend(group.split(','))
        self.assertEqual(set(keys), set(entries))
        first = list(dict.fromkeys(keys))
        self.assertEqual(first, entries)

    def test_referenced_local_paths_exist(self):
        for relative in re.findall(r'\\path\{([^}]+)\}', self.text):
            resolved = (ROOT/relative).exists() if '/' in relative else any(ROOT.rglob(relative))
            self.assertTrue(resolved, relative)

    def test_no_external_tex_inputs(self):
        self.assertIsNone(re.search(r'\\(?:input|include|externaldocument)\s*\{', self.text))
        self.assertEqual([p.name for p in ROOT.glob('*.tex')], ['supplementary_material_gamma.tex'])
        self.assertEqual([p.name for p in ROOT.glob('*.pdf')], ['supplementary_material_gamma.pdf'])

    def test_package_version_consistency(self):
        version = (ROOT/'VERSION').read_text().strip()
        self.assertRegex(version, r'^\d+\.\d+\.\d+$')
        citation = (ROOT/'CITATION.cff').read_text()
        match = re.search(r'^version:\s*(\S+)\s*$', citation, re.M)
        self.assertIsNotNone(match)
        self.assertEqual(match.group(1), version)
        zenodo = json.loads((ROOT/'.zenodo.json').read_text())
        self.assertEqual(zenodo['version'], version)
        self.assertEqual(len(zenodo['creators']), 5)
        self.assertEqual(zenodo['language'], 'eng')

    def test_lotfi_certificate(self):
        data = json.loads((ROOT/'results'/'lotfi_orders.json').read_text())
        methods = data['fields']['coupled_quadratic']['methods']
        self.assertEqual([methods[n]['leading_degree'] for n in ('M6','M9','M12')], [5,7,9])
        self.assertIn(r'\label{supp:prop:lotfi-orders}', self.text)


if __name__ == '__main__':
    unittest.main()
