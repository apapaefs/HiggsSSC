"""VBF_H contracts with synthetic LHE files; no event generation."""
from dataclasses import replace
import gzip
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import run_pipeline
from hfourlepton import higher_order as ho, workflow
from hfourlepton import run_four_lepton_campaign as four
from hfourlepton import run_powheg_vbf as vbf
from test_fourlepton_ho import cfg


def fixture(ebeam=6500., masses=(125.0011342, 124.98, 125.05)):
    card = vbf.input_card(ebeam)
    parameters, weights = card.split('<initrwgt>', 1)
    events = []
    for index, mass in enumerate(masses):
        weight = (6., -1., 4.)[index % 3]
        event = (f'<event>\n5 9999 {weight} 80.398 .007 .118\n'
                 '2 -1 0 0 501 0 0 0 100 100 0 0 9\n'
                 '1 -1 0 0 502 0 0 0 -100 100 0 0 9\n'
                 f'25 1 1 2 0 0 0 0 0 {mass} {mass} 0 9\n'
                 '2 1 1 2 501 0 25 0 30 39.051248 0 0 9\n'
                 '1 1 1 2 502 0 -25 0 -30 39.051248 0 0 9\n')
        event += '<rwgt>\n' + ''.join(f"<wgt id='{i}'> {weight} </wgt>\n" for i in range(1001, 1008))
        events.append(event + '</rwgt>\n</event>\n')
    return ('<LesHouchesEvents version="3.0">\n<!--\n' + '\n'.join(vbf.PROVENANCE) + '\n' +
            parameters + '-->\n<header>\n<initrwgt>' + weights + '</header>\n<init>\n' +
            f'2212 2212 {ebeam} {ebeam} -1 -1 -1 -1 -4 1\n-1 .1 6 9999\n</init>\n' +
            ''.join(events) + '</LesHouchesEvents>\n')


class VBFTests(unittest.TestCase):
    def test_both_energies_native_rate_disjoint_decays_and_finite_width(self):
        for energy in (13, 40):
            with self.subTest(energy=energy), tempfile.TemporaryDirectory() as tmp, patch.object(four, 'SCRIPT_DIR', Path(tmp)):
                source = Path(tmp) / 'vbf.lhe'
                source.write_text(fixture(energy * 500, masses=tuple(125. + i * .001 for i in range(9))))
                original = source.read_bytes()
                settings = replace(cfg('--ebeam', str(energy * 500), '--nevents', '3',
                    '--vbf-lhe', str(source)), dry_run=False)
                samples = [s for s in ho.all_samples() if s.production_mode == 'vbf']
                info = ho.inspect_production(source, samples[0], settings)
                self.assertEqual(info['cross_section_pb'], 3.)
                self.assertEqual(info['negative_events'], 3)
                slices = []
                for sample in samples:
                    sliced = ho.obtain_production(sample, samples, settings)
                    slices.append(gzip.open(sliced, 'rt').read().split('<event>')[1])
                    self.assertEqual(four.validate_stable_higgs_lhe(sliced, settings, 3, sample=sample), (125., .004152))
                    self.assertEqual(four.parse_lhe_init_cross_section_pb(sliced), 3.)
                    decayed = Path(tmp) / 'decayed.lhe'
                    decayed.write_text(source.read_text().replace('-1 .1 6 9999', '.001 .0001 6 9999'))
                    normalized = ho.normalization(sample, decayed, settings, sliced)
                    self.assertEqual(normalized['normalization_kind'], 'native_lhe')
                    self.assertAlmostEqual(normalized['target_cross_section_pb'], 3 * four.signal_target_branching_fraction(sample, settings))
                    self.assertNotIn('ihixs', normalized)
                self.assertEqual(len(set(slices)), 3)
                self.assertEqual(source.read_bytes(), original)

    def test_invalid_source_provenance_beams_pdf_width_scale_weights_fail(self):
        changes = [(vbf.PROVENANCE[0], 'some Hjj process'), (vbf.COMMIT, 'wrong revision'),
            ('hwidth 0.004152', 'hwidth 0.005'), ('higgsfixedwidth 1', 'higgsfixedwidth 0'),
            ('renscfact 1', 'renscfact 2'), ('lhans1 335900', 'lhans1 336100'),
            ('2212 2212 6500.0 6500.0', '2212 2212 20000 20000'),
            ('125.0011342', '125.2'), ('80.398 .007 .118', 'nan .007 .118'),
            ("<wgt id='1001'> 6.0 </wgt>", ''),
            ("<wgt id='1001'> 6.0 </wgt>", "<wgt id='1001'> nan </wgt>")]
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'bad.lhe'
            for old, new in changes:
                with self.subTest(change=old):
                    self.assertIn(old, fixture())
                    source.write_text(fixture().replace(old, new))
                    with self.assertRaises(ValueError):
                        vbf.inspect_lhe(source, SimpleNamespace(ebeam=6500., higgs_mass=125.), 3)

    def test_vbf_cannot_accidentally_enter_mcatnlo_or_ihixs_path(self):
        sample = next(s for s in ho.all_samples() if s.production_mode == 'vbf')
        with self.assertRaisesRegex(ValueError, 'VBF uses POWHEG'):
            ho.nlo_cards(sample, cfg(), Path('/unused'), 3)
        card = ho.herwig_card(sample, Path('/synthetic/vbf.lhe'), 1, 3, cfg())
        self.assertIn('MaxPtIsMuF Yes', card)
        self.assertNotIn('InitialStateReconOption Rapidity', card)
        self.assertIn('NNPDF40_nlo_as_01180_qed', card)
        self.assertEqual(four.sample_manifest_entry(sample, cfg())['hard_model'], 'SM')
        self.assertEqual(four.sample_manifest_entry(ho.all_samples()[0], cfg())['hard_model'], 'pure_heft')

    def test_vbf_setup_is_idempotent_and_checks_its_component(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            obj = vbf.shared.Pipeline.__new__(vbf.shared.Pipeline)
            obj.work, obj.powheg = root, root / 'POWHEG-BOX-V2'
            obj.config = {'cc': '/usr/bin/gcc', 'cxx': '/usr/bin/g++', 'fc': '/usr/bin/gfortran'}
            process = obj.powheg / 'VBF_H'
            calls = []
            def call(label, command, cwd=None):
                calls.append((label, command))
                if label == 'powheg-vbf-source':
                    process.mkdir(parents=True)
                    (process / '.git').write_text('synthetic submodule')
                elif label == 'powheg-vbf-build':
                    (process / 'pwhg_main').write_bytes(b'synthetic binary')
            obj.call, obj.capture = call, lambda command: vbf.COMMIT
            with patch.object(vbf, 'source_identity', return_value={'source': 'synthetic'}):
                obj.setup_powheg_vbf()
                obj.setup_powheg_vbf()
                self.assertEqual([x[0] for x in calls], ['powheg-vbf-source', 'powheg-vbf-build'])
                self.assertEqual((process / 'svn.version').read_text(), vbf.version_include())
                self.assertIn('F77=/usr/bin/gfortran -fno-automatic -Jobj-gfortran -Iobj-gfortran', calls[1][1])
                obj.record_component('powheg_vbf', {'source': 'changed'})
                with self.assertRaisesRegex(ValueError, 'installation differs'):
                    obj.setup_powheg_vbf()

    def test_preparation_includes_vbf_card_and_independent_vh(self):
        for energy in (13, 40):
            with self.subTest(energy=energy), tempfile.TemporaryDirectory() as tmp:
                args = run_pipeline.parse_args(['prepare', '--workflow', 'fourlepton-ho',
                    '--energy-tev', str(energy), '--run-tag', 'vbf_cards'])
                with patch.object(workflow, 'REPO', Path(tmp)), patch.object(workflow, 'configuration', return_value={}):
                    plan = workflow.prepare(args)
                self.assertEqual((plan / 'signal_vbf_h_stable/powheg.input').read_text(), vbf.input_card(energy * 500))
                self.assertFalse((plan / 'signal_vbf_h_stable/process.mg5').exists())
                self.assertTrue((plan / 'signal_z_h_stable/compile.mg5').exists())
                inventory = json.loads((plan / 'plan.json').read_text())
                self.assertEqual(len(inventory['samples']), 27)

    def test_resume_rejects_changed_configuration_and_completed_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            arguments = ['--powheg-dir', str(root / 'source'), '--run-dir', str(root / 'run'),
                         '--nevents', '3', '--ebeam', '6500', '--seed', '12']
            def produce(command, *args, **kwargs):
                (root / 'run' / vbf.LHE_NAME).write_text(fixture())
            with patch.object(vbf, 'source_identity', return_value={'synthetic': True}), \
                 patch.object(vbf.mg5_guard, 'run_guarded', side_effect=produce) as run, \
                 patch.object(vbf.driver, 'log_is_complete', return_value=True):
                vbf.main(arguments)
                vbf.main(arguments)
                self.assertEqual(run.call_count, 1)
                with self.assertRaisesRegex(ValueError, 'configuration changed'):
                    vbf.main(arguments + ['--nevents', '6'])
                (root / 'run' / vbf.LHE_NAME).write_text('corrupted')
                with self.assertRaisesRegex(ValueError, 'completed VBF LHE changed'):
                    vbf.main(arguments)


if __name__ == '__main__':
    unittest.main()
