import copy
import gzip
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from hgammagamma import run_loop_gammagamma as loop
from hgammagamma import run_gammagamma_ho_campaign as ho
from hgammagamma import evaluate_loop_gammagamma as evaluate
from hgammagamma import frozen_gammagamma_classifier as frozen
import analyze_lo_varfiles as analysis


def lhe_text(process=loop.PROCESS, incoming=(21,21), outgoing=(22,22), weights=(1.,1.,1.)):
    records=[]
    for weight in weights:
        rows=[f"4 1 {weight} 100 .007 .118"]
        for pdg in incoming:
            rows.append(f"{pdg} -1 0 0 501 502 0 0 100 100 0 0 9")
        for i,pdg in enumerate(outgoing):
            rows.append(f"{pdg} 1 1 2 0 0 {30 if i==0 else -30} 0 0 30 0 0 9")
        records.append("<event>\n"+"\n".join(rows)+"\n</event>\n")
    return (f'<LesHouchesEvents version="3.0">\n<header><MG5ProcCard>\ngenerate {process}\n</MG5ProcCard></header>\n'
            '<init>\n2212 2212 20000 20000 0 0 335900 335900 3 1\n5.0 .1 1.0 1\n</init>\n'
            + ''.join(records)+'</LesHouchesEvents>\n')


class LoopGammaGammaTest(unittest.TestCase):
    def args(self,path,*extra):
        return ho.parse_args(['--output-dir',str(path),'--nevents','3','--no-herwig-module',*extra])

    def test_adapter_restores_nominal_samples_and_dispatch(self):
        original=ho.SAMPLES
        inspector=ho.inspect_lhe
        with loop.loop_adapter():
            self.assertEqual(ho.SAMPLES,(loop.SAMPLE,))
            self.assertIs(ho.inspect_lhe,loop.inspect_lhe)
        self.assertIs(ho.SAMPLES,original)
        self.assertIs(ho.inspect_lhe,inspector)
        self.assertEqual(len(ho.SAMPLES),4)

    def test_process_excludes_higgs_and_shower_has_no_mcatnlo_or_signal_decay(self):
        args=self.args('/tmp/loop-test','--mg5-fortran','/usr/bin/gfortran','--mg5-cxx','/usr/bin/g++')
        process,settings=loop.mg5_cards(args)
        self.assertIn('generate g g > a a / h [noborn=QCD]',process)
        self.assertIn('output_dependencies internal --no_save',process)
        self.assertIn('set fortran_compiler /usr/bin/gfortran --no_save',process)
        self.assertIn('set cpp_compiler /usr/bin/g++ --no_save',process)
        self.assertNotIn('output /tmp/loop-test -f',process)
        card=loop.herwig_card(args)
        self.assertNotIn('LongTransBoost',card)
        self.assertNotIn('SelectDecayModes',card)
        self.assertNotIn('BranchingRatio',card)
        self.assertIn('AllowedToReOpen No',card)
        self.assertIn('NNPDF40_nlo_as_01180_qed',card)
        values=json.loads(settings)
        self.assertEqual(values['lhaid'],335900)
        self.assertEqual(values['pta'],10.)
        self.assertEqual(values['etaa'],6.)
        self.assertEqual(values['dynamical_scale_choice'],3)

    def test_prepared_manifest_is_explicitly_lo_and_unity_weight(self):
        with TemporaryDirectory() as tmp,loop.loop_adapter():
            args=self.args(tmp)
            manifest=loop.prepare(args)
            self.assertEqual(manifest['matching'],'LO')
            self.assertEqual(manifest['hard_accuracy'],'loop-induced LO + PS')
            self.assertEqual(manifest['weight_scale'],1.)
            self.assertTrue(manifest['physics_scope']['higgs_excluded'])
            self.assertEqual(manifest,json.loads((ho.sample_dir(args,loop.SAMPLE)/'campaign.json').read_text()))

    def test_native_loop_lhe_needs_no_mcatnlo_provenance(self):
        with TemporaryDirectory() as tmp:
            path=Path(tmp)/'events.lhe.gz'
            with gzip.open(path,'wt') as f: f.write(lhe_text())
            result=loop.inspect_lhe(path,loop.SAMPLE,self.args(tmp))
            self.assertEqual(result['events'],3)
            self.assertEqual(result['cross_section_pb'],5.)
            self.assertTrue(result['process_validation']['higgs_excluded'])
            self.assertEqual(result['negative_events'],0)
            with self.assertRaisesRegex(ValueError,'MC@NLO'):
                ho.inspect_lhe(path,ho.SAMPLES[1],self.args(tmp))

    def test_process_topology_and_weight_mismatches_fail(self):
        cases=[dict(process='g g > a a [noborn=QCD]'),dict(incoming=(1,-1)),
               dict(outgoing=(22,25)),dict(weights=(1.,-1.,1.))]
        for change in cases:
            with self.subTest(change=change),TemporaryDirectory() as tmp:
                path=Path(tmp)/'events.lhe'
                path.write_text(lhe_text(**change))
                with self.assertRaises(ValueError): loop.inspect_lhe(path,loop.SAMPLE,self.args(tmp))

    def test_run_card_patch_controls_hidden_pdf_and_scale_defaults(self):
        with TemporaryDirectory() as tmp:
            path=Path(tmp)/'run_card.dat'
            path.write_text(' 230000 = lhaid ! legacy PDF\n 20 = pta\n')
            settings=loop.run_settings(self.args(tmp))
            loop.patch_run_card(path,settings)
            text=path.read_text()
            self.assertIn('335900 = lhaid',text)
            self.assertIn('lhapdf = pdlabel1',text)
            self.assertIn('3 = dynamical_scale_choice',text)
            self.assertIn('10.0 = pta',text)
            loop.patch_run_card(path,settings)
            self.assertEqual(path.read_text(),text)

    def test_native_mass_interval_is_parsed_written_and_audited(self):
        with TemporaryDirectory() as tmp:
            args=loop.parse_args(['--output-dir',tmp,'--nevents','3','--no-herwig-module',
                '--gen-mgg-min','50','--gen-mgg-max','70'])
            self.assertEqual((loop.run_settings(args)['mmaa'],loop.run_settings(args)['mmaamax']),(50.,70.))
            path=Path(tmp)/'events.lhe';path.write_text(lhe_text())
            result=loop.inspect_lhe(path,loop.SAMPLE,args)
            self.assertEqual(result['born_mgg_validation']['minimum_observed_gev'],60.)
            args.gen_mgg_min=61.
            with self.assertRaisesRegex(ValueError,'mass interval'):loop.inspect_lhe(path,loop.SAMPLE,args)
        with self.assertRaisesRegex(ValueError,'mass interval'):
            loop.parse_args(['--gen-mgg-min','100','--gen-mgg-max','90'])

    def test_compiled_export_copy_has_independent_writes_and_rebased_internal_links(self):
        with TemporaryDirectory() as tmp:
            source=Path(tmp).resolve()/'smoke/mg5_process';destination=Path(tmp).resolve()/'target/mg5_process'
            (source/'Cards').mkdir(parents=True);(source/'bin/internal').mkdir(parents=True)
            (source/'Cards/proc_card_mg5.dat').write_text(f'import model {loop.MODEL}\ngenerate {loop.PROCESS}\n')
            for name in ['Cards/run_card.dat','Cards/param_card.dat','bin/internal/banner.py','bin/internal/madevent_interface.py']:
                (source/name).write_text('original\n')
            (source.parent/'campaign.json').write_text(json.dumps({'sample':loop.SAMPLE.name,
                'process':loop.PROCESS,'completed':['analyze'],'shower_completion':{},
                'generator_artifacts':{name:loop.four.file_record(source/'Cards'/name)
                    for name in ('run_card.dat','param_card.dat')}}))
            (source/'absolute_link').symlink_to(source/'Cards/param_card.dat')
            with patch.object(loop.completion,'validate_completion') as validate:
                evidence=loop.clone_completed_export(source,destination)
            validate.assert_called_once()
            self.assertIn('no hard links',evidence['copy_policy'])
            self.assertEqual((destination/'absolute_link').resolve(),destination/'Cards/param_card.dat')
            (destination/'Cards/run_card.dat').write_text('new cuts\n')
            self.assertEqual((source/'Cards/run_card.dat').read_text(),'original\n')
            (source/'Cards/param_card.dat').write_text('changed mass\n')
            with patch.object(loop.completion,'validate_completion'), self.assertRaisesRegex(ValueError,'cards changed'):
                loop.clone_completed_export(source,Path(tmp)/'another/mg5_process')

    def test_loop_accuracy_requires_completion_in_varfile_loader(self):
        with TemporaryDirectory() as tmp:
            directory=Path(tmp)/loop.SAMPLE.name
            directory.mkdir()
            (directory/'campaign.json').write_text(json.dumps({
                'hard_accuracy':'loop-induced LO + PS','run_tag':'pilot',
                'sample':loop.SAMPLE.name,'completed':['generate']}))
            with self.assertRaisesRegex(ValueError,'matching completed analysis'):
                analysis.ho_sample_metadata(directory,'pilot',directory/'summary.dat',{},None,'Backgrounds')

    def test_mc_error_counts_failed_showers_as_zero_and_groups_response_rows(self):
        # Four attempted hard events; the first contributes two correlated
        # response rows of weight 1/2, the second fails the cuts, two fail showering.
        error=evaluate.source_mc_error([True,True,False],[.5,.5,1.],[0,0,1],
            attempted=4,cross_section=8.,denominator=4.)
        self.assertEqual(error,2.)

    def test_positive_source_bound_does_not_treat_zero_events_as_zero_background(self):
        result=evaluate.positive_source_upper_bound([False,False,False],[.4,.6,1.],[0,0,1],
            attempted=100,cross_section=10.)
        self.assertEqual(result['selected_positive_sources'],0)
        self.assertAlmostEqual(result['selected_cross_section_upper_pb'],10.*(1.-.05**.01))
        result=evaluate.positive_source_upper_bound([True,True,False],[.4,.6,1.],[0,0,1],
            attempted=100,cross_section=10.)
        self.assertEqual(result['selected_positive_sources'],1)
        self.assertGreater(result['selected_cross_section_upper_pb'],.1)
        from scipy.stats import binomtest
        exact=binomtest(1,100,alternative='less').proportion_ci(confidence_level=.95,method='exact')
        self.assertAlmostEqual(result['selected_cross_section_upper_pb'],10.*exact.high)
        with self.assertRaisesRegex(ValueError,'unit total response'):
            evaluate.positive_source_upper_bound([True,True],[.7,.7],[0,0],attempted=1,cross_section=10.)

    def test_frozen_cut_summary_rejects_repeat_loop_and_classifier_summary(self):
        record={'metadata':{'mode':'cuts','detector_response':'ssc','luminosity_fb':100.,
            'cuts':[{'variable':'m_gg','min':120.,'max':130.}]},'samples':[
                {'sample':'bkg_prompt_aa','category':'Backgrounds','selected_cross_section_pb':1.,
                 'expected_events':100000.,'detector_response':'ssc'},
                {'sample':'signal_gg_h_aa','category':'Signal','selected_cross_section_pb':.1,
                 'expected_events':10000.,'detector_response':'ssc'}]}
        with TemporaryDirectory() as tmp:
            path=Path(tmp)/'summary.json'
            path.write_text(json.dumps(record))
            self.assertEqual(evaluate.frozen_nominal(path),record)
            invalid=copy.deepcopy(record)
            invalid['samples'].append({**invalid['samples'][0],'sample':loop.SAMPLE.name})
            path.write_text(json.dumps(invalid))
            with self.assertRaisesRegex(ValueError,'already included'): evaluate.frozen_nominal(path)
            invalid=copy.deepcopy(record)
            invalid['metadata']['mode']='xgboost'
            path.write_text(json.dumps(invalid))
            with self.assertRaisesRegex(ValueError,'SSC cut summary'): evaluate.frozen_nominal(path)

    def test_frozen_classifier_preserves_feature_order_and_candidate_only_preselection(self):
        import numpy as np
        class Model:
            def predict_proba(self, values):
                self.values=values
                return np.asarray([[.001,.999]])
        model=Model()
        eligible_row={name:float(i+1) for i,name in enumerate(frozen.FEATURE_NAMES)}
        eligible_row.update(m_gg=180.,n_selected_photons=2.)
        missing={name:-999. for name in frozen.FEATURE_NAMES}
        missing['n_selected_photons']=1.
        eligible,scores=frozen.score_rows(model,[eligible_row,missing])
        self.assertEqual(eligible.tolist(),[True,False])
        self.assertEqual(model.values.tolist(),[[eligible_row[name] for name in frozen.FEATURE_NAMES]])
        self.assertEqual(scores[0],.999)
        self.assertTrue(np.isnan(scores[1]))
        self.assertTrue((eligible & (scores>=.992))[0])

    def test_frozen_classifier_rejects_changed_feature_order(self):
        metadata=dict(mode='xgboost',status='complete',detector_response='ssc',
            evaluation_partition='test',feature_names=frozen.FEATURE_NAMES,best_threshold=.992)
        frozen.validate_metadata(metadata)
        with self.assertRaisesRegex(ValueError,'feature order'):
            frozen.validate_metadata({**metadata,'feature_names':list(reversed(frozen.FEATURE_NAMES))})


if __name__=='__main__': unittest.main()
