import gzip
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from hfourlepton import run_reducible_controls as controls


def lhe(count=5, weight=1.0):
    header = '<LesHouchesEvents version="3.0">\n<init>\n2212 2212 20000 20000 0 0 331900 331900 3 1\n2.0 0.01 1.0 1\n</init>\n'
    events = ''.join(f'<event>\n1 1 {weight} 125 0.007 0.118\n11 1 0 0 0 0 {i} 0 0 {i+1} 0 0 9\n</event>\n' for i in range(count))
    return header + events + '</LesHouchesEvents>\n'


class ReducibleControlsTest(unittest.TestCase):
    def test_partition_covers_all_records_and_preserves_native_rate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root/'source.lhe'
            source.write_text(lhe())
            result = controls.split_lhe(source, root/'parts', 2)
            self.assertEqual([(s['start'], s['stop']) for s in result['shards']], [(0,2),(2,4),(4,5)])
            self.assertEqual(result['source_events'], 5)
            original = [data for kind,data in controls.lhe_parts(source) if kind=='event']
            rebuilt=[]
            for shard in result['shards']:
                rebuilt.extend(data for kind,data in controls.lhe_parts(shard['lhe']) if kind=='event')
                self.assertEqual(controls.four.parse_lhe_init_cross_section_pb(Path(shard['lhe'])), 2.)
            self.assertEqual(rebuilt, original)
            self.assertEqual(controls.split_lhe(source, root/'parts', 2), result)
            source.write_text(lhe(6))
            with self.assertRaisesRegex(ValueError, 'changed'):
                controls.split_lhe(source, root/'parts', 2)

    def test_partition_rejects_modified_artifact_and_incomplete_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); source = root/'source.lhe'; source.write_text(lhe())
            result=controls.split_lhe(source,root/'parts',2)
            Path(result['shards'][0]['lhe']).write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError,'changed'):
                controls.split_lhe(source,root/'parts',2)
            source.write_text(lhe().replace('</LesHouchesEvents>',''))
            with self.assertRaisesRegex(ValueError,'incomplete'):
                controls.split_lhe(source,root/'other',2)

    def test_unit_source_weights_required_for_pooling(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); source=root/'source.lhe'; source.write_text(lhe().replace('1 1 1.0 125','1 1 -1.0 125',1))
            with self.assertRaisesRegex(ValueError,'positive unit'):
                controls.split_lhe(source,root/'parts',2)

    def test_finite_card_preserves_decay_and_places_eof_guard_before_save(self):
        original='set A B\n# forced decay configuration\nsaverun sample theGenerator\n'
        result=controls.finite_card(original)
        self.assertIn('# forced decay configuration',result)
        self.assertLess(result.index('AllowedToReOpen No'),result.index('saverun'))
        self.assertEqual(result.count('set /Herwig/Generators/theGenerator:DumpPeriod -1'),1)
        self.assertLess(result.index('DumpPeriod -1'),result.index('saverun'))
        with self.assertRaises(ValueError): controls.finite_card(result)

    def test_actual_cuts_are_read_from_export(self):
        with tempfile.TemporaryDirectory() as tmp:
            card=Path(tmp)/'run_card.dat'
            card.write_text(' 20.0 = ptb ! b jet pT\n -1 = etab\n 0.4 = drbb\n 7 = unrelated\n')
            self.assertEqual(controls.actual_generation_cuts(card),{'ptb':'20.0','etab':'-1','drbb':'0.4'})

    def test_inclusive_policy_sets_hidden_jet_cuts_and_preserves_model_inputs(self):
        original='20.0 = ptj ! jet cut\n 0.4 = drjl\n 331900 = lhaid\n 5 = maxjetflavor\n'
        modified=controls.rewrite_inclusive_run_card(original)
        self.assertIn('331900 = lhaid',modified)
        self.assertIn('5 = maxjetflavor',modified)
        self.assertEqual(controls.rewrite_inclusive_run_card(modified),modified)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); (root/'Cards').mkdir(); (root/'Source').mkdir()
            (root/'Cards/run_card.dat').write_text(modified)
            compiled=''.join(f'      {key.upper()} = {value:.8f}D+00\n' for key,value in controls.INCLUSIVE_HEAVY_FLAVOUR_CUTS.items())
            (root/'Source/run_card.inc').write_text(compiled+'      MAXJETFLAVOR = 5\n')
            result=controls.audit_inclusive_generation(root)
            self.assertEqual(result['effective_generation_cuts']['maxjetflavor'],5.)
            (root/'Source/run_card.inc').write_text(compiled.replace('PTJ = 0.00000000D+00','PTJ = 20.00000000D+00'))
            with self.assertRaisesRegex(ValueError,'ptj'):
                controls.audit_inclusive_generation(root)
        with self.assertRaisesRegex(ValueError,'duplicate'):
            controls.rewrite_inclusive_run_card('20 = ptj\n10 = ptj\n')

    def test_generation_policy_is_scoped_and_changes_production_identity(self):
        def identity(*args,**kwargs): return {'base':'unchanged'}
        with mock.patch.object(controls.four,'production_stage_identity',identity):
            before=controls.four.patch_run_card
            with controls.apply_generation_policy(False):
                self.assertEqual(controls.four.production_stage_identity(),{'base':'unchanged'})
                self.assertIs(controls.four.patch_run_card,before)
            with controls.apply_generation_policy(True):
                sample=controls.four.SAMPLES_BY_NAME['bkg_zbb_ee']
                self.assertEqual(controls.four.production_stage_identity(sample)['control_generation_policy'],controls.generation_policy(True,sample.name))
                self.assertIsNot(controls.four.patch_run_card,before)
            self.assertIs(controls.four.patch_run_card,before)
            self.assertIs(controls.four.production_stage_identity,identity)

    def test_flavour_creation_changes_only_zbb_incoming_definition(self):
        original='import model sm\ndefine p = '+controls.STANDARD_INCOMING+'\ndefine j = '+controls.STANDARD_INCOMING+'\ngenerate p p > e+ e- b b~\n'
        modified=controls.flavour_creation_process_card(original,'bkg_zbb_ee')
        self.assertIn('define p = '+controls.LIGHT_INCOMING+'\n',modified)
        self.assertIn('define j = '+controls.STANDARD_INCOMING+'\n',modified)
        self.assertEqual(controls.flavour_creation_process_card(original,'bkg_ttbar'),original)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); card=root/'process.dat'; source=root/'source.lhe'
            card.write_text(modified)
            text=lhe(1).replace('1 1 1.0 125','3 1 1.0 125').replace('11 1 0 0 0 0 0 0 0 1 0 0 9',
                '21 -1 0 0 501 0 0 0 100 100 0 0 9\n21 -1 0 0 0 501 0 0 -100 100 0 0 9\n11 1 1 2 0 0 0 0 0 1 0 0 9')
            source.write_text(text)
            audit=controls.audit_incoming_partons(source,card,'bkg_zbb_ee')
            self.assertEqual(audit['incoming_parton_audit']['pairs'],{'21,21':1})
            source.write_text(text.replace('21 -1','5 -1',1))
            with self.assertRaisesRegex(ValueError,'incoming bottom'):
                controls.audit_incoming_partons(source,card,'bkg_zbb_ee')

    def test_native_lo_card_passes_strict_source_reader_audit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); source=root/'source.lhe'; source.write_text(lhe())
            template=Path('hfourlepton/HW-template.in').read_text()
            card=root/'card.in'
            card.write_text(controls.finite_card(template.replace('LHEFILE.lhe.gz',str(source))))
            self.assertEqual(controls.completion._audit_source_card(card),str(source.resolve()))

    def test_branching_evidence_preserves_runtime_precision(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); log=root/'read.log'; card=root/'card.in'; evidence=root/'evidence.json'
            log.write_text('Herwig runtime\n0.106832\n0.106049\n0.106832\n0.106049\n')
            card.write_text('runtime card')
            with mock.patch.object(controls.completion,'prompt_branching_transform',create=True,return_value={'proof':True}):
                self.assertEqual(controls.prompt_evidence(log,card,evidence),{'proof':True})
            data=json.loads(evidence.read_text())
            self.assertEqual(list(data['branching_ratio_tokens'].values()).count('0.106832'),2)
            self.assertAlmostEqual(data['selected_branching_sums']['24'],0.212881)
            log.write_text('0.1\n')
            with self.assertRaisesRegex(ValueError,'exactly four'):
                controls.prompt_evidence(log,card,evidence)

    def test_branching_evidence_integrates_with_completion_transform(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); log=root/'read.log'; card=root/'card.in'; evidence=root/'evidence.json'
            log.write_text('0.106832\n0.106049\n0.106832\n0.106049\n')
            card.write_text('\n'.join(controls.four.PROMPT_DECAY_COMMANDS['w_emu'])+'\n'
                'create HiggsSSC::LHEBranchingRatioReweighter /Herwig/Generators/PromptBRReweighter\n'
                'insert X:PostHadronizationHandlers 0 /Herwig/Generators/PromptBRReweighter\n'
                +''.join(f'get {mode}:BranchingRatio\n' for mode in controls.PROMPT_MODES))
            transform=controls.prompt_evidence(log,card,evidence)
            self.assertAlmostEqual(transform['factor'],0.212881**2)
            self.assertLess(transform['factor_bounds'][0],transform['factor'])
            self.assertGreater(transform['factor_bounds'][1],transform['factor'])


if __name__=='__main__': unittest.main()
