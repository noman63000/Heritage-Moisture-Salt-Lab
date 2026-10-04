"""Regression tests for the explicit 1.2.1 overnight QA -> spin-up -> freeze workflow."""
from __future__ import annotations
import tempfile, unittest
from pathlib import Path
from unittest.mock import patch
from hmsl.workflow import prepare_selected_cases
from hmsl.common import InputError

ROOT=Path(__file__).resolve().parents[1]

class DummyProject:
    def __init__(self,path):
        self.path=Path(path)
        self.spin_ok=False
    def readiness(self):
        return {'spinup_valid':self.spin_ok,
                'spinup':{'baseline_run':'R00'} if self.spin_ok else {}}

class Overnight121Tests(unittest.TestCase):
    def test_automatic_queue_runs_qa_then_spinup_then_freeze_and_continues_failure(self):
        with tempfile.TemporaryDirectory() as td:
            p=DummyProject(td);events=[]
            def qa(project,rid,hours,progress,cancel):
                events.append(('qa',rid))
                if rid=='R04': raise InputError('deliberate QA fail')
                return {'passed':True}
            def sp(project,rid,progress,cancel):
                events.append(('spinup',rid));project.spin_ok=True;return {'passed':True}
            def fr(project,qa_reference_run=None):
                events.append(('freeze',qa_reference_run));return {'passed':True}
            with patch('hmsl.workflow.require_trial',lambda project:None), \
                 patch('hmsl.workflow.numerical_qa',qa), \
                 patch('hmsl.workflow.spinup',sp), \
                 patch('hmsl.workflow.valid_case_qa',lambda project,rid:rid!='R04'), \
                 patch('hmsl.workflow.freeze',fr):
                result=prepare_selected_cases(p,['R03','R04','R05'],baseline_run='R00')
            self.assertEqual(events,[('qa','R03'),('qa','R04'),('qa','R05'),('spinup','R00'),('freeze','R03')])
            self.assertEqual(result['qa_passed'],2)
            self.assertEqual(result['qa_failed'],1)
            self.assertEqual(result['spinup'],'PASS')
            self.assertEqual(result['freeze'],'PASS')

    def test_interface_explicitly_says_one_button_runs_all_three_stages(self):
        text=(ROOT/'web/index.html').read_text(encoding='utf-8')
        self.assertIn('Overnight preparation: QA → Spin-up → Freeze',text)
        self.assertIn('One button performs all three stages automatically.',text)
        self.assertIn("Run QA → Spin-up → Freeze overnight",text)
        self.assertIn('You do not need to press the manual Spin-up or Freeze buttons',text)
        self.assertIn('Stage 2/3 Water/heat spin-up',text)
        self.assertIn('Stage 3/3 Freeze reduced inputs',text)

if __name__=='__main__': unittest.main(verbosity=2)
