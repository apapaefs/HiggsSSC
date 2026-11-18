import array
import csv
import math
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import numpy as np
from hfourlepton.reducible_isolation_sidebands import second_mass
from hgammagamma.dijet_response import moments
from hgammagamma.evaluate_dijet_minbias import generated_cross_section
from hgammagamma.run_dijet_bound import set_cuts, check_include


class DiagnosticsTests(unittest.TestCase):
    def test_minbias_cross_section_units_and_last_digit_error(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "minbias.out"
            path.write_text("(nb)\nTotal (from generated events): 10000 10000 90.2(5)e+06\n"
                            "Events carry unit weights.\n")
            self.assertEqual(generated_cross_section(path), (10000, 90.2e6, 0.5e6))
            path.write_text(path.read_text().replace("(nb)", "(pb)"))
            with self.assertRaises(ValueError):
                generated_cross_section(path)

    def test_scored_pairs_group_before_squaring_and_cap_per_event(self):
        dtype = [("source", int), ("eligible_jets", int), ("f1", float), ("f2", float),
                 ("fq", float), ("fg", float)]
        pairs = np.array([(0, 3, .1, .2, .6, .1), (0, 3, .3, .2, .6, .1),
                          (1, 2, .2, .2, .6, .1)], dtype=dtype)
        mask = np.ones(3, dtype=bool)
        sw, sw2, n = moments(pairs, mask, 3, 1, "independent_pair_sum")
        self.assertAlmostEqual(sw, .12)
        self.assertAlmostEqual(sw2, .08**2 + .04**2)
        self.assertEqual(n, 2)
        self.assertEqual(moments(pairs, mask, 3, 1, "all_quark_any_correlation"),
                         (1.6, 1.0 + .6**2, 2))
        self.assertEqual(moments(pairs, ~mask, 3, 1, "all_quark_independent"), (0., 0., 0))
        self.assertAlmostEqual(moments(pairs, mask, 3, 1, "all_gluon_two_eligible")[0], .01)

    def test_loose_quartet_pairing_nearest_z_and_unpaired(self):
        a = dict(lepton_pt=np.array([[45.,45.,8.,8.]]*2), lepton_eta=np.zeros((2,4)),
                 lepton_phi=np.array([[0.,math.pi,0.,math.pi]]*2),
                 lepton_charge=np.array([[-1,1,-1,1],[1,1,-1,1]]),
                 lepton_flavour=np.array([[11,11,13,13]]*2))
        values = second_mass(a)
        self.assertAlmostEqual(values[0],16.)
        self.assertTrue(np.isnan(values[1]))

    def test_effective_bin_cuts_reject_different_upper_edge(self):
        with tempfile.TemporaryDirectory() as tmp:
            card=Path(tmp)/"card"
            card.write_text(" 20.0 = ptj ! cut\n")
            set_cuts(card,{"ptj":5.,"ptjmax":15.})
            self.assertIn("5.0 = ptj",card.read_text())
            self.assertIn("15.0 = ptjmax",card.read_text())
            card.write_text(" PTJMAX = 30d0\n")
            with self.assertRaises(ValueError):
                check_include(card,{"ptjmax":15.})

    @unittest.skipUnless(shutil.which("root-config"),"requires ROOT")
    def test_dijet_source_moments_and_pair_envelopes(self):
        import ROOT
        root=Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory() as tmp:
            tmp=Path(tmp)
            flags=subprocess.check_output(["root-config","--cflags","--libs"],text=True).split()
            analyzer=tmp/"analyzer"
            subprocess.run(["c++","-std=c++17",root/"hgammagamma/LOAnalysis/Code/HwSimPostAnalysis_gammagamma_SSC.cc",
                            *flags,"-o",analyzer],check=True)
            path=tmp/"events.root"
            file=ROOT.TFile(str(path),"RECREATE")
            tree=ROOT.TTree("Data","test SSC jet pairs")
            leaves={}
            for name in ("numJets","numbJets","numoutgoing","theJets","thebJets","partons","evweight"):
                if name.startswith("num"):
                    leaf=array.array("i",[0]); fmt="/I"
                elif name=="evweight":
                    leaf=array.array("d",[1.]); fmt="/D"
                else:
                    leaf=array.array("d",[0.]*500); fmt="[5][100]/D"
                tree.Branch(name,leaf,name+fmt)
                leaves[name]=leaf
            for i,(eta,phi,flavour) in enumerate(((.5,0.,21),(-.5,math.pi,1),(-.5,math.pi-.45,1))):
                pt=62.5/math.cosh(eta)
                for j,value in enumerate((62.5,pt*math.cos(phi),pt*math.sin(phi),pt*math.sinh(eta),flavour)):
                    leaves["theJets"][100*j+i]=value
                    leaves["partons"][100*j+i]=value
            for count in (2,3):
                leaves["numJets"][0]=count
                leaves["numoutgoing"][0]=count
                tree.Fill()
            tree.Write()
            file.Close()
            out=tmp/"result.csv"
            subprocess.run([analyzer,"--dijet-bound",path,out,"23"],check=True)
            with out.open() as f:
                r={(v["scenario"],float(v["fake_scale"])):v for v in csv.DictReader(f)}
            def val(k,scale=1): return float(r[k,scale]["sum_probability"])
            self.assertEqual(int(r["two_eligible_only",1]["contributing_events"]),1)
            self.assertEqual(int(r["independent_pair_sum",1]["contributing_events"]),2)
            self.assertLess(val("two_eligible_only"),val("independent_pair_sum"))
            self.assertLess(val("independent_pair_sum"),val("all_quark_independent"))
            self.assertLess(val("all_quark_independent"),val("all_quark_any_correlation"))
            self.assertAlmostEqual(val("all_quark_independent",2)/val("all_quark_independent"),4.)
