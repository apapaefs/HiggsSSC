"""Synthetic LHE encoding/decay-boundary checks; no generators or showers."""
import gzip
import re
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest

from hfourlepton import run_four_lepton_campaign as four


def fixture(*, decayed=False, scaled=False):
    header = ('<LesHouchesEvents version="3.0">\n<header><initrwgt>\n'
              '<weightgroup name="scale"><weight id="20">zeta</weight>\n'
              '<weight id="3">alpha</weight></weightgroup>\n</initrwgt></header>\n'
              '<init>\n2212 2212 6500 6500 0 0 336100 336100 -4 1\n'
              f'{1 if decayed else 10} 0 1 1\n</init>\n')
    events = []
    for index, weight in enumerate((10., -2., 0.)):
        factor = .1 if decayed else 1.
        optional_factor = factor if scaled else 1.
        particle = f'25 {2 if decayed else 1} 0 0 0 0 {index} 0 0 125 125 0 9\n'
        daughters = '11 1 1 1 0 0 0 0 0 0 0 0 9\n' if decayed else ''
        events.append(f'<event>\n{2 if decayed else 1} 1 {weight * factor} 125 .007 .118\n'
                      + particle + daughters
                      + f'<weights> {2 * weight * optional_factor} {-weight * optional_factor} </weights>\n</event>\n')
    return header + ''.join(events) + '</LesHouchesEvents>\n'


def write_pair(root, *, scaled=False):
    source = root / 'source.lhe'
    decayed = root / 'decayed.lhe.gz'
    source.write_text(fixture())
    with gzip.open(decayed, 'wt') as stream:
        stream.write(fixture(decayed=True, scaled=scaled))
    return source, decayed


class FourLeptonWeightTests(unittest.TestCase):
    def test_compact_values_follow_declaration_order_before_lexical_hwsim_mapping(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / 'source.lhe'
            path.write_text(fixture())
            definitions = four.parse_lhe_optional_weight_definitions(path)
            payload = next(payload for event, payload in four.lhe_blocks(path) if event)
            weights, compact = four.event_optional_weights(payload, definitions, str(path))
            self.assertTrue(compact)
            self.assertEqual(weights, {'20': 20., '3': -10.})
            mapping = four.lhe_optional_weight_mapping(path, definitions)
            self.assertEqual([item['id'] for item in mapping['mapping']], ['3', '20'])
            self.assertEqual(mapping['compact_events'], 3)

    def test_invalid_compact_payloads_fail_closed(self):
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / 'source.lhe'
            path.write_text(fixture())
            definitions = four.parse_lhe_optional_weight_definitions(path)
            for bad in ('<weights>1</weights>', '<weights>1 2 3</weights>',
                        '<weights>1 2', '<weights/>',
                        '<weights>1 2</weights><weights>1 2</weights>',
                        '<weights>1 2</weights><rwgt><wgt id="20">1</wgt></rwgt>',
                        '<weights>nan 2</weights>', '<weights>inf 2</weights>',
                        '<weights>word 2</weights>'):
                with self.subTest(bad=bad):
                    path.write_text('<event>\n1 1 1 1 1 1\n' + bad + '\n</event>\n')
                    with self.assertRaises(ValueError):
                        four.lhe_optional_weight_mapping(path, definitions)
            path.write_text('<event>\n1 1 1 1 1 1\n<weights>1 2</weights>\n</event>\n')
            with self.assertRaisesRegex(ValueError, 'count differs'):
                four.lhe_optional_weight_mapping(path, [])
            path.write_text('<event>\n1 1 1 1 1 1\n<weights>nan 2</weights>\n</event>\n')
            diagnostic = four.lhe_optional_weight_mapping(path, definitions, allow_nonfinite=True)
            self.assertEqual(diagnostic['nonfinite_event_counts'], {'20': 1})

    def test_compact_madspin_copy_scales_variations_and_preserves_every_other_byte(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            source, decayed = write_pair(root)
            before = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in (source, decayed)}
            cfg, sample = SimpleNamespace(run_dir=root / 'run'), SimpleNamespace(name='signal')
            output, proof = four.prepare_shower_lhe_weights(decayed, source, sample, cfg)
            self.assertEqual(proof['identity']['decay_factor'], .1)
            self.assertEqual(proof['rescaled_optional_values'], 4)
            with gzip.open(output, 'rt') as stream:
                prepared = stream.read()
            original = fixture(decayed=True)
            self.assertEqual(re.sub(r'<rwgt>.*?</rwgt>', '', prepared, flags=re.S),
                             re.sub(r'<weights>.*?</weights>', '', original, flags=re.S))
            definitions = four.parse_lhe_optional_weight_definitions(output)
            values = [four.event_optional_weights(payload, definitions, str(output))[0]
                      for event, payload in four.lhe_blocks(output) if event]
            self.assertEqual(values, [{'20': 2., '3': -1.}, {'20': -.4, '3': .2}, {'20': 0., '3': 0.}])
            self.assertEqual(four.lhe_optional_weight_mapping(output, definitions)['compact_events'], 0)
            for path, expected in before.items():
                self.assertEqual((path.read_bytes(), path.stat().st_mtime_ns), expected)
            output_before = (output.read_bytes(), output.stat().st_mtime_ns)
            self.assertEqual(four.prepare_shower_lhe_weights(decayed, source, sample, cfg)[0], output)
            self.assertEqual((output.read_bytes(), output.stat().st_mtime_ns), output_before)
            source.write_text(source.read_text().replace('20.0 -10.0', '21.0 -10.0'))
            with self.assertRaisesRegex(ValueError, 'differ from saved inputs'):
                four.prepare_shower_lhe_weights(decayed, source, sample, cfg)

    def test_already_scaled_compact_values_are_not_scaled_twice(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            source, decayed = write_pair(root, scaled=True)
            output, proof = four.prepare_shower_lhe_weights(decayed, source,
                SimpleNamespace(name='signal'), SimpleNamespace(run_dir=root / 'run'))
            self.assertEqual(proof['rescaled_optional_values'], 0)
            definitions = four.parse_lhe_optional_weight_definitions(output)
            payload = next(payload for event, payload in four.lhe_blocks(output) if event)
            self.assertEqual(four.event_optional_weights(payload, definitions, str(output))[0], {'20': 2., '3': -1.})

    def test_background_conversion_needs_no_decay_factor(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            source, _ = write_pair(root)
            output, proof = four.prepare_shower_lhe_weights(source, None,
                SimpleNamespace(name='background'), SimpleNamespace(run_dir=root / 'run'))
            self.assertEqual(proof['identity']['decay_factor'], 1.)
            self.assertEqual(proof['rescaled_optional_values'], 0)
            definitions = four.parse_lhe_optional_weight_definitions(output)
            payload = next(payload for event, payload in four.lhe_blocks(output) if event)
            self.assertEqual(four.event_optional_weights(payload, definitions, str(output))[0], {'20': 20., '3': -10.})

    def test_unexplained_weights_changed_momenta_and_event_count_are_rejected(self):
        for bad, expected in (
            (fixture(decayed=True).replace('20.0 -10.0', '17.0 -10.0'), 'unexplained'),
            (fixture(decayed=True).replace('20.0 -10.0', '2.0 -10.0'), 'unexplained'),
            (fixture(decayed=True).replace('2 1 1.0 125', '2 1 2.0 125'), 'nominal weights'),
            (fixture(decayed=True).replace('25 2 0 0 0 0 0 0 0', '25 2 0 0 0 0 5 0 0'), 'momenta'),
            (fixture(decayed=True).replace('25 2 0 0 0 0', '24 2 0 0 0 0'), 'identity'),
            (re.sub(r'<event>.*?</event>\n', '', fixture(decayed=True), count=1, flags=re.S), 'nominal weights'),
            (re.sub(r'<event>\n2 1 0.0.*?</event>\n', '', fixture(decayed=True), flags=re.S), 'event counts'),
        ):
            with self.subTest(expected=expected), TemporaryDirectory() as tmp:
                root = Path(tmp)
                source, decayed = write_pair(root)
                with gzip.open(decayed, 'wt') as stream:
                    stream.write(bad)
                with self.assertRaisesRegex(ValueError, expected):
                    four.prepare_shower_lhe_weights(decayed, source,
                        SimpleNamespace(name='signal'), SimpleNamespace(run_dir=root / 'run'))
                self.assertFalse((root / 'run/weight_inputs/signal/events.lhe.gz').exists())


if __name__ == '__main__':
    unittest.main()
