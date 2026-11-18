import copy
import json
from pathlib import Path
import tempfile
import unittest

from hfourlepton import hf_decay_audit as audit
from hfourlepton import run_weighted_reducible_pilot as pilot
from hfourlepton import run_reducible_controls as controls
from hgammagamma import ho_shower_completion as completion
from tests.test_reducible_controls import lhe


class WeightedPilotTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name).resolve()
        self.table = self.root / 'decays.in'
        self.table.write_text(''.join(
            f'decaymode {parent}->D0,nu_e,e+; 0.2 1 Decayer\n'
            f'decaymode {parent}->pi+,pi-; 0.8 1 Decayer\n'
            for parent in controls.four.HEAVY_MESON_PARENTS))
        self.log = self.root / 'run.log'

    def row(self, number, semileptonic):
        factor = .4 if semileptonic else 1.6
        return {'event': number, 'before': 1., 'factor': factor, 'after': factor,
                'parents': [{'pdg': 511, 'mode': '/Herwig/Particles/B0/B0->D0,nu_e,e+;' if semileptonic else '/Herwig/Particles/B0/B0->pi+,pi-;',
                             'original_sum': 1., 'proposal_sum': 1.6, 'mode_bias': 4. if semileptonic else 1., 'factor': factor}]}

    def write_rows(self, rows):
        self.log.write_text('\n'.join(audit.PREFIX + json.dumps(row) for row in rows))

    def test_independent_mode_weight_and_nonselected_tail(self):
        self.write_rows([self.row(1, True), self.row(2, False)])
        summary = audit.inspect_log(self.log, self.table)
        self.assertEqual(summary['events'], 2)
        self.assertAlmostEqual(summary['root_moments']['sum_weight'], 2.)
        self.assertAlmostEqual(summary['root_moments']['sum_weight_squared'], 2.72)
        self.assertEqual(summary['parent_counts']['B0'], {'all': 2, 'selected': 1})

    def test_changed_mode_or_factor_or_replayed_event_fails(self):
        for field, value in [('factor', 1.), ('proposal_sum', 1.), ('mode', '/Herwig/Particles/B0/B0->*MatchAny;')]:
            row = self.row(1, True)
            row['parents'][0][field] = value
            self.write_rows([row])
            with self.assertRaises(ValueError):
                audit.inspect_log(self.log, self.table)
        self.write_rows([self.row(1, True), self.row(1, True)])
        with self.assertRaisesRegex(ValueError, 'repeated'):
            audit.inspect_log(self.log, self.table)

    def test_partonic_colour_connection_canonicalization(self):
        self.table.write_text(self.table.read_text().replace('B0->pi+,pi-; 0.8', 'B0->d=dbar,dbar=d; 8.0D-1'))
        row = self.row(1, False)
        row['parents'][0]['mode'] = '/Herwig/Particles/B0/B0->d=dbar=dbar,d;'
        self.write_rows([row])
        summary = audit.inspect_log(self.log, self.table)
        self.assertEqual(summary['parent_counts']['B0'], {'all':1,'selected':0})

    def test_final_weight_cannot_hide_wrong_product(self):
        row = self.row(1, True)
        row['parents'].append(copy.deepcopy(row['parents'][0]))
        self.write_rows([row])
        with self.assertRaisesRegex(ValueError, 'product'):
            audit.inspect_log(self.log, self.table)

    def test_root_closure_and_zero_response_source_denominator(self):
        self.write_rows([self.row(1, True), self.row(3, False)])
        plugin, source, card = (self.root/name for name in ('plugin.so','plugin.cc','card.in'))
        plugin.write_text('private compiled plugin')
        source.write_text('private audited source')
        commands, _, _ = controls.four.semileptonic_decay_commands(self.table)
        card.write_text('\n'.join([f'read {self.table}', f'globallibrary {plugin}',
                                   'create HiggsSSC::LHEHeavyFlavorBranchingRatioReweighter handler',
                                   'insert generator:PostDecayHandlers 0 handler', *commands]))
        evidence = audit.build_audit(self.log, self.table, card, source, plugin)
        record = {'sample':'bkg_zbb_ee_hfbiased','attempted_events':3,'saved_events':2,
                  'heavy_flavour_weight_audit':evidence,
                  'artifacts':{'herwig_card':completion._artifact(card), 'run_log':completion._artifact(self.log)}}
        source_weights = {key:3. for key in completion.MOMENTS}
        source_weights['negative_events'] = 0
        root_weights = evidence['summary']['root_moments']
        audit.validate_audit(record, source_weights, root_weights, verify_file=True)
        with self.assertRaisesRegex(ValueError,'every saved'):
            audit.validate_audit({**record,'saved_events':1},source_weights,root_weights)
        with self.assertRaisesRegex(ValueError,'ROOT moments'):
            audit.validate_audit(record,source_weights,{**root_weights,'sum_weight':3.})
        self.table.write_text(self.table.read_text().replace('0.2','0.3'))
        with self.assertRaisesRegex(ValueError,'changed'):
            audit.validate_audit(record,source_weights,root_weights,verify_file=True)

    def test_exact_disjoint_partial_parent_slices(self):
        source = self.root/'parent.lhe'
        source.write_text(lhe(10))
        partition = controls.split_lhe(source, self.root/'parent_parts', 4)
        record = {'source_events': 10, 'shards': partition['shards']}
        slices = pilot.select_range(record, 3, 5, self.root/'pilot_parts')
        self.assertEqual([(row['start'],row['stop']) for row in slices], [(3,4),(4,8)])
        original = [data for kind,data in controls.lhe_parts(source) if kind=='event']
        actual = [data for row in slices for kind,data in controls.lhe_parts(row['lhe']) if kind=='event']
        self.assertEqual(actual, original[3:8])
        for row in slices:
            self.assertEqual(controls.four.parse_lhe_init_cross_section_pb(Path(row['lhe'])), 2.)
        with self.assertRaisesRegex(ValueError, 'outside'):
            pilot.select_range(record, 8, 3, self.root/'invalid')

    def test_changed_parent_shard_is_rejected(self):
        source = self.root/'parent.lhe'; source.write_text(lhe(4))
        partition = controls.split_lhe(source, self.root/'parts', 4)
        partition['shards'][0]['sha256'] = '0'*64
        with self.assertRaisesRegex(ValueError, 'changed'):
            pilot.select_range({'source_events':4,'shards':partition['shards']},0,4,self.root/'new')


if __name__ == '__main__':
    unittest.main()
