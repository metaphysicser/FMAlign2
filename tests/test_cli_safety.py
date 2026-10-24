"""Synthetic regressions for input errors, argument errors and workspace ownership."""
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

binary = Path(sys.argv.pop(1)).resolve()
mock = Path(__file__).with_name('mock_aligner.py').resolve()


class CliSafety(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='fmalign2-safety-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.work = self.root / 'work space'
        (self.work / 'nested').mkdir(parents=True)
        self.preserved = {'keep.txt': b'preserve this', 'tiny.fasta': b'existing input',
                          'aligned.fasta': b'existing output', 'nested/keep.txt': b'nested data'}
        for name, content in self.preserved.items():
            (self.work / name).write_bytes(content)
        self.prefix = 'call ' if os.name == 'nt' else ''
        self.template = self.prefix + f'"{sys.executable}" "{mock}" "{{input}}" "{{output}}"'

    def run_case(self, content, options=None, template=None):
        (self.root / 'input.fasta').write_text(content)
        (self.root / 'command.txt').write_text(self.template if template is None else template)
        args = {'-i': 'input.fasta', '-o': 'output.fasta', '-p': 'command.txt',
                '-w': self.work.name, '-t': '1', '-l': '30'}
        args.update(options or {})
        result = subprocess.run([str(binary)] + [v for pair in args.items() for v in pair],
                                cwd=self.root, capture_output=True, text=True, timeout=30)
        for name, original in self.preserved.items():
            self.assertTrue((self.work / name).is_file(), (name, result.stderr))
            self.assertEqual((self.work / name).read_bytes(), original)
        self.assertEqual(list(self.work.glob('probe-*')), [], 'Probe directory leaked')
        self.assertNotIn('terminate called', result.stderr)
        return result

    def reject(self, content, options=None, template=None, message='Error:'):
        result = self.run_case(content, options, template)
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn(message, result.stderr)
        self.assertFalse((self.root / 'output.fasta').exists())

    def test_existing_directory_success(self):
        result = self.run_case('>a\nACGT\n>b\nACGA\n')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue((self.root / 'output.fasta').exists())

    def test_backend_failure(self):
        self.reject('>a\nACGT\n', template=self.prefix +
                    f'"{sys.executable}" -c "import sys; sys.exit(7)" {{input}} {{output}}')

    def test_backend_missing_output(self):
        self.reject('>a\nACGT\n', template=self.prefix +
                    f'"{sys.executable}" -c "import sys; sys.exit(0)" {{input}} {{output}}')

    def test_missing_placeholder_cleanup(self):
        self.reject('>a\nACGT\n', template='echo {input}')

    def test_fastq_errors(self):
        good = '@a\nACGT\n+\n!!!!\n'
        for bad in ('@b\nACGA\n+\n!!\n', '@b\nACGA\n+\n!!!!!\n', '@b\nACGA\n+\n'):
            with self.subTest(bad=bad):
                self.reject(good + bad, message='Invalid FASTA/FASTQ input')

    def test_fastq_success(self):
        result = self.run_case('@a\nACGT\n+\n!!!!\n@b\nACGA\n+\n!!!!\n')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((self.root / 'output.fasta').read_text(), '>a\nACGT\n>b\nACGA\n')

    def test_invalid_parameters(self):
        invalid = {'-f': ['invalid'], '-v': ['2', '-1', '1junk'],
                   '-t': ['0', '-1', 'abc', '2junk', '9999999999999999999999'],
                   '-l': ['0', '-1', 'abc', '30junk', '9999999999999999999999']}
        for option, values in invalid.items():
            for value in values:
                with self.subTest(option=option, value=value):
                    self.reject('>a\nACGT\n', {option: value}, message='Invalid ' + option)

    def test_discovered_anchors(self):
        rng = random.Random(20261024)
        sequence = ''.join(rng.choice('ACGT') for _ in range(160))
        for mode in ('accurate', 'fast'):
            with self.subTest(mode=mode):
                result = self.run_case(f'>a\n{sequence}\n>b\n{sequence}\n', {'-f': mode})
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertNotIn('Sequence divide parts: 1 ', result.stdout)
                self.assertEqual((self.root / 'output.fasta').read_text(),
                                 f'>a\n{sequence}\n>b\n{sequence}\n')

    def test_missing_parameter_value(self):
        self.reject('>a\nACGT\n', {'-t': ''}, message='Missing value for argument: -t')


if __name__ == '__main__':
    unittest.main()
