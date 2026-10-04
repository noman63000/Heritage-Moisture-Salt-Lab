"""Workflow/UI regression tests for 1.1.1 queue behaviour."""
from __future__ import annotations
import json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from hmsl.workflow import prepare_selected_cases
from hmsl.common import InputError
ROOT=Path(__file__).resolve().parents[1]

class DummyProject:
    def __init__(self,path):self.path=Path(path)
    def readiness(self):return {'spinup_valid':True,'spinup':{'baseline_run':'R00'}}

class QueueTests(unittest.TestCase):
    def test_prepare_queue_continues_after_one_case_fails(self):
        with tempfile.TemporaryDirectory() as td:
            p=DummyProject(td);called=[]
            def qa(project,rid,hours,progress,cancel):
                called.append(rid)
                if rid=='R04':raise InputError('deliberate QA failure')
                return {'passed':True}
            frozen=[]
            with patch('hmsl.workflow.require_trial',lambda project:None), \
                 patch('hmsl.workflow.numerical_qa',qa), \
                 patch('hmsl.workflow.valid_case_qa',lambda project,rid:rid!='R04'), \
                 patch('hmsl.workflow.freeze',lambda project,qa_reference_run=None:frozen.append(qa_reference_run) or {'passed':True}):
                result=prepare_selected_cases(p,['R03','R04','R05'])
            self.assertEqual(called,['R03','R04','R05'])
            self.assertEqual(result['qa_passed'],2);self.assertEqual(result['qa_failed'],1)
            self.assertEqual(result['freeze'],'PASS');self.assertEqual(frozen,['R03'])
            saved=json.loads((Path(td)/'prep_queue.json').read_text())
            self.assertEqual(saved['case_qa'][1]['run_id'],'R04');self.assertEqual(saved['case_qa'][1]['qa'],'FAIL')

    def test_dashboard_has_persistent_multi_selection_and_range(self):
        text=(ROOT/'web/index.html').read_text(encoding='utf-8')
        self.assertIn('localStorage.setItem(selectionKey()',text)
        self.assertIn("job('prep_queue')",text)
        self.assertIn('Select range',text)
        self.assertIn('One failed case does not stop the remaining queue',text)

if __name__=='__main__':unittest.main(verbosity=2)
