"""Integration contract tests; no bioinformatics dependencies or input data needed.

Backends are substituted only at the subprocess boundary. These tests verify
runHiC orchestration, not numerical correctness of the supplied algorithms.
"""
import contextlib
import io
import logging
import os
from pathlib import Path
import runpy
import subprocess
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from runHiC import binning as backend


class FastBinningTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='runhic space ')
        self.root = Path(self.temp.name) / '矩阵测试'
        self.root.mkdir()
        self.pairs_dir = self.root / 'filtered-test'
        self.pairs_dir.mkdir()
        self.pairs = self.pairs_dir / 'sample-filtered.pairs.gz'
        self.pairs.write_bytes(b'pairs fixture')
        Path(str(self.pairs) + '.px2').write_bytes(b'index fixture')
        (self.pairs_dir / 'sample.completed').touch()
        self.base = self.root / 'sample.10.cool'
        self.mcool = self.root / 'sample.mcool'
        self.utilities = types.ModuleType('runHiC.utilities')

        def chromsizes(_):
            sizes = self.root / 'genome sizes.txt'
            sizes.write_text('chr1\t100\n')
            return str(sizes), 'test'

        self.utilities.chromsizes_from_pairs = chromsizes
        self.utilities.cleanFile = lambda p: Path(p).unlink(missing_ok=True)
        self.modules = patch.dict(sys.modules, {'runHiC.utilities': self.utilities})
        self.modules.start()
        self.cli = runpy.run_path(str(ROOT / 'scripts' / 'runHiC'), run_name='runhic_test_cli')

    def tearDown(self):
        self.modules.stop()
        self.temp.cleanup()

    def parse(self, *args):
        with patch.object(sys, 'argv', ['runHiC', *args]):
            return self.cli['getargs']()

    @staticmethod
    def fake_backend(command):
        if command[1].endswith('fast_pairix_to_cool.py'):
            Path(command[-1]).write_bytes(b'base matrix fixture')
        else:
            Path(command[command.index('--out') + 1]).write_bytes(b'mcool fixture')

    def build(self, **kwargs):
        backend.mcool_from_pairs(str(self.pairs), str(self.base), str(self.mcool), '10,50,100', **kwargs)

    def test_resolution_contract(self):
        self.assertEqual(backend.resolve_resolutions('100,10,50,10'), (10, '10,50,100'))
        self.assertEqual(backend.resolve_resolutions([5000, 10000]), (5000, '5000,10000'))
        self.assertEqual(backend.resolve_resolutions('50,100', 10), (10, '10,50,100'))
        self.assertEqual(backend.resolve_resolutions('b', 10), (10, 'B'))
        self.assertEqual(backend.resolve_resolutions('N', 10), (10, 'N'))
        for spec in ('B', '0,10', '-10,20', '10,25', '10,,50', 'x', ''):
            with self.subTest(spec=spec), self.assertRaises(ValueError):
                backend.resolve_resolutions(spec)

    def test_shared_parser_defaults_and_pileup_options(self):
        binargs, _ = self.parse('binning', '-f', str(self.pairs_dir))
        pileargs, _ = self.parse('pileup', '--nproc', '3', '--ice-threads', '2',
                                 '--resolutions', 'B', '--base-resolution', '10',
                                 '--ignore-diags', '0', '--min-nnz', '4')
        self.assertEqual(binargs.resolutions, backend.DEFAULT_RESOLUTIONS)
        self.assertEqual((pileargs.nproc, pileargs.ice_threads, pileargs.threads), (3, 2, 8))
        self.assertEqual((pileargs.ignore_diags, pileargs.min_nnz), (0, 4))

    def test_invalid_options_fail_before_work(self):
        options = [('--resolutions', 'B'), ('--resolutions', '10,25'),
                   ('--memory-gb', 'nan'), ('--tol', 'inf'), ('--nproc', '0'),
                   ('--min-nnz', '-1'), ('--compression', '10'), ('--window-bp', '1')]
        for option in options:
            with self.subTest(option=option), contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as exc:
                self.parse('binning', '-f', str(self.pairs_dir), *option)
            self.assertEqual(exc.exception.code, 2)

    def test_backend_arguments_paths_and_cleanup(self):
        with patch.object(backend.subprocess, 'check_call', side_effect=self.fake_backend) as call:
            self.build(nproc=3, ice_threads=2, window_bp=5000000, buffer_contacts=4000000,
                       temp_dir=str(self.root / 'scratch space'), ignore_diags=0, ignore_dist=30,
                       min_nnz=3, min_count=4, mad_max=6, blacklist=str(self.root / 'mask file.bed'),
                       tol=1e-6, max_iters=500, memory_gb=160, chunksize=101, io_chunk=17,
                       compression=0, convergence_policy='store_final')
        loader, zoom = [c.args[0] for c in call.call_args_list]
        self.assertTrue(all(isinstance(c.args[0], list) and not c.kwargs for c in call.call_args_list))
        self.assertEqual(loader[0], sys.executable)
        self.assertEqual(loader[-2], str(self.pairs))
        self.assertTrue(loader[-3].endswith(':10'))
        for command, options in ((loader, {'--nproc':'3', '--window-bp':'5000000', '--buffer-contacts':'4000000', '--io-chunk':'17', '--compression':'0'}),
                                 (zoom, {'--threads':'2', '--resolutions':'10,50,100', '--ignore-diags':'0', '--ignore-dist':'30', '--min-nnz':'3', '--min-count':'4', '--mad-max':'6', '--tol':'1e-06', '--max-iters':'500', '--memory-gb':'160', '--chunksize':'101', '--convergence-policy':'store_final'})):
            for key, value in options.items():
                self.assertEqual(command[command.index(key) + 1], value)
            self.assertEqual(command[command.index('--temp-dir') + 1], str(self.root / 'scratch space'))
        self.assertEqual(zoom[zoom.index('--blacklist') + 1], str(self.root / 'mask file.bed'))
        self.assertFalse(self.base.exists())
        self.assertFalse(Path(str(self.base) + '.runhic.json').exists())
        self.assertTrue(self.mcool.exists())

    def test_nonconvergence_retains_base_and_resume_skips_loading(self):
        def fail_zoom(command):
            self.fake_backend(command)
            if command[1].endswith('fast_zoomify.py'):
                raise subprocess.CalledProcessError(2, command)
        with patch.object(backend.subprocess, 'check_call', side_effect=fail_zoom), self.assertRaises(subprocess.CalledProcessError) as exc:
            self.build()
        self.assertEqual(exc.exception.returncode, 2)
        self.assertTrue(self.base.exists())
        self.assertTrue(Path(str(self.base) + '.runhic.json').exists())
        with patch.object(backend.subprocess, 'check_call', side_effect=self.fake_backend) as call:
            self.build(resume=True, max_iters=500, keep_base_cool=True)
        self.assertEqual(call.call_count, 1)
        self.assertIn('--resume', call.call_args.args[0])
        self.assertTrue(self.base.exists())

    def test_resume_after_failure_before_mcool_creation(self):
        def fail_zoom(command):
            if command[1].endswith('fast_zoomify.py'):
                raise subprocess.CalledProcessError(1, command)
            self.fake_backend(command)
        with patch.object(backend.subprocess, 'check_call', side_effect=fail_zoom), self.assertRaises(subprocess.CalledProcessError):
            self.build()
        self.assertFalse(self.mcool.exists())
        with patch.object(backend.subprocess, 'check_call', side_effect=self.fake_backend) as call:
            self.build(resume=True)
        self.assertEqual(call.call_count, 1)

    def test_changed_input_blocks_resume(self):
        with patch.object(backend.subprocess, 'check_call', side_effect=self.fake_backend):
            self.build(keep_base_cool=True)
        self.pairs.write_bytes(b'changed pairs input')
        with patch.object(backend.subprocess, 'check_call') as call, self.assertRaisesRegex(ValueError, 'changed'):
            self.build(resume=True)
        call.assert_not_called()

    def test_missing_index_and_existing_output_fail_before_loading(self):
        index = Path(str(self.pairs) + '.px2')
        index.unlink()
        with patch.object(backend.subprocess, 'check_call') as call, self.assertRaises(FileNotFoundError):
            self.build()
        call.assert_not_called()
        index.touch()
        self.mcool.write_bytes(b'previous result')
        with patch.object(backend.subprocess, 'check_call') as call, self.assertRaises(FileExistsError):
            self.build()
        call.assert_not_called()
        self.assertEqual(self.mcool.read_bytes(), b'previous result')

    def test_deprecated_option_warns_without_forwarding(self):
        with patch.object(backend.subprocess, 'check_call', side_effect=self.fake_backend) as call, self.assertLogs(backend.log, level='WARNING'):
            self.build(max_split=2, nproc=3)
        for item in call.call_args_list:
            self.assertNotIn('--max-split', item.args[0])
        zoom = call.call_args.args[0]
        self.assertEqual(zoom[zoom.index('--threads') + 1], '3')

    def test_pipeline_marker_and_lock_on_success(self):
        args, commands = self.parse('binning', '-f', str(self.pairs_dir), '--resolutions', '10,50,100')
        previous = os.getcwd()
        try:
            os.chdir(self.root)
            with patch.object(backend.subprocess, 'check_call', side_effect=self.fake_backend):
                self.cli['binning'](args, commands)
        finally:
            os.chdir(previous)
        output = self.root / 'coolers-test'
        self.assertTrue((output / 'sample-filtered.completed').exists())
        self.assertTrue((output / 'sample-filtered.mcool').exists())
        self.assertFalse((output / 'sample-filtered.lock').exists())
        self.assertFalse((output / 'sample-filtered.10.cool').exists())

    def test_pipeline_failure_cleans_lock_without_marker(self):
        args, commands = self.parse('binning', '-f', str(self.pairs_dir), '--resolutions', '10')
        previous = os.getcwd()
        try:
            os.chdir(self.root)
            with patch.object(backend, 'mcool_from_pairs', side_effect=subprocess.CalledProcessError(2, ['zoom'])), self.assertRaises(subprocess.CalledProcessError):
                self.cli['binning'](args, commands)
        finally:
            os.chdir(previous)
        output = self.root / 'coolers-test'
        self.assertFalse((output / 'sample-filtered.completed').exists())
        self.assertFalse((output / 'sample-filtered.lock').exists())

    def test_pileup_preserves_user_matrix_settings(self):
        args, commands = self.parse('pileup', '--nproc', '3', '--min-nnz', '4', '--ignore-diags', '0')
        captured = []
        def record(current, _):
            captured.append((current.nproc, current.min_nnz, current.ignore_diags))
        namespace = self.cli['pileup'].__globals__
        with patch.dict(namespace, mapping=record, filtering=record, binning=record):
            self.cli['pileup'](args, commands)
        self.assertEqual(captured, [(3, 4, 0)] * 3)

    def test_loader_failure_stops_before_zoom_and_provenance(self):
        with patch.object(backend.subprocess, 'check_call', side_effect=subprocess.CalledProcessError(1, ['loader'])) as call, self.assertRaises(subprocess.CalledProcessError):
            self.build()
        self.assertEqual(call.call_count, 1)
        self.assertFalse(self.mcool.exists())
        self.assertFalse(Path(str(self.base) + '.runhic.json').exists())

    def test_cli_preserves_nonconvergence_exit_status(self):
        argv = ['runHiC', 'binning', '-f', str(self.pairs_dir), '--resolutions', '10',
                '--logFile', str(self.root / 'failure.log')]
        previous = os.getcwd()
        logger = logging.getLogger()
        old_handlers, old_level = list(logger.handlers), logger.level
        try:
            os.chdir(self.root)
            with patch.object(sys, 'argv', argv), patch.object(backend, 'mcool_from_pairs', side_effect=subprocess.CalledProcessError(2, ['zoom'])), self.assertRaises(SystemExit) as exc:
                runpy.run_path(str(ROOT / 'scripts' / 'runHiC'), run_name='__main__')
            self.assertEqual(exc.exception.code, 2)
        finally:
            os.chdir(previous)
            for handler in list(logger.handlers):
                if handler not in old_handlers:
                    logger.removeHandler(handler)
                    handler.close()
            logger.setLevel(old_level)


if __name__ == '__main__':
    unittest.main()
