"""Repeatable implementation checks. These are not site validation tests."""
from __future__ import annotations
import copy, json, shutil, tempfile, unittest
from pathlib import Path
import numpy as np
import pandas as pd
from hmsl.common import InputError, Cancelled, atomic_json, sha256
from hmsl.climate import read_climate, ordinal_hours, Climate
from hmsl.material import read_material, Material, parse_m6
from hmsl.project import Project, new_project
from hmsl.solver import ColumnModel, State, Solubility
from hmsl.workflow import run_transport, freeze, numerical_qa, spinup
ROOT=Path(__file__).resolve().parents[1]

class Sandbox(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.cfg=json.loads((ROOT/'projects/demo/project.json').read_text())
        self.obj=json.loads((ROOT/'projects/demo/materials/SYNTHETIC.json').read_text())
    def tearDown(self):self.tmp.cleanup()
    def climate(self,dates=None,hours=None,**updates):
        n=len(dates) if dates is not None else 24
        df=pd.DataFrame({'date':dates or ['2020-01-01']*n,'hour_utc':hours if hours is not None else list(range(n)),
          'tas_C':25.,'hurs_pct':60.,'pr_mm_h':0.,'rsds_W_m2':0.,'rlds_W_m2':380.,'sfcWind_m_s':2.})
        for k,v in updates.items():df[k]=v
        p=self.root/'climate.csv';df.to_csv(p,index=False);return p
    def demo(self,name):
        dst=self.root/name
        shutil.copytree(ROOT/'projects/demo',dst,ignore=shutil.ignore_patterns('outputs','*.npz','qa.json','spinup.json','last_job.json','last_error.log'))
        return Project(dst)

class InputTests(Sandbox):
    def test_gregorian_leap_day(self):
        p=self.climate(['2020-02-28','2020-02-29'],[23,0]);self.assertEqual(read_climate(p).n,2)
    def test_noleap_skips_gregorian_leap_day(self):
        p=self.climate(['2040-02-28','2040-03-01'],[23,0]);self.assertEqual(read_climate(p,'365_day').n,2)
        with self.assertRaises(InputError):read_climate(p,'gregorian')
    def test_noleap_rejects_feb29(self):
        p=self.climate(['2040-02-28','2040-02-29'],[23,0])
        with self.assertRaises(InputError):read_climate(p,'365_day')
    def test_360_day_feb30(self):
        p=self.climate(['2040-02-30','2040-03-01'],[23,0]);self.assertEqual(read_climate(p,'360_day').n,2)
    def test_duplicate_hour_rejected(self):
        p=self.climate(['2020-01-01']*2,[0,0])
        with self.assertRaises(InputError):read_climate(p)
    def test_gap_rejected(self):
        p=self.climate(['2020-01-01']*2,[0,2])
        with self.assertRaises(InputError):read_climate(p)
    def test_nan_hour_rejected(self):
        p=self.climate(['2020-01-01']*2,[0,np.nan])
        with self.assertRaises(InputError):read_climate(p)
    def test_kelvin_not_silently_converted(self):
        with self.assertRaises(InputError):read_climate(self.climate(tas_C=298.15))
    def test_RH_invalid_rejected(self):
        with self.assertRaises(InputError):read_climate(self.climate(hurs_pct=101))
    def test_fractional_RH_flagged(self):
        self.assertTrue(read_climate(self.climate(hurs_pct=.6)).warnings)
    def test_negative_rain_rejected(self):
        with self.assertRaises(InputError):read_climate(self.climate(pr_mm_h=-.1))
    def test_nan_numeric_rejected(self):
        with self.assertRaises(InputError):read_climate(self.climate(rlds_W_m2=np.nan))
    def test_checksum_mismatch_rejected(self):
        with self.assertRaises(InputError):read_climate(self.climate(),expected_sha256='wrong')
    def test_audit_workbook_not_mistaken_for_forcing(self):
        p=self.root/'audit.xlsx';pd.DataFrame({'note':['not climate']}).to_excel(p,sheet_name='Audit',index=False)
        with self.assertRaises(InputError):read_climate(p)
    def test_Climate_sheet_accepted(self):
        csv=self.climate();p=self.root/'climate.xlsx';pd.read_csv(csv).to_excel(p,sheet_name='Climate',index=False)
        self.assertEqual(read_climate(p).n,24)

class MaterialTests(Sandbox):
    def test_missing_source_stops(self):
        with self.assertRaisesRegex(InputError,'Missing material'):read_material(self.root/'absent.m6')
    def test_m6_units_and_log_function(self):
        p=self.root/'SPECIFICATION_FIXTURE.m6'
        p.write_text('[STORAGE_BASE_PARAMETERS]\nRHO = 1800 kg/m3\nCE = 850 J/kgK\nTHETA_POR = 0.38 m3/m3\nTHETA_EFF = 0.34 m3/m3\n[TRANSPORT_BASE_PARAMETERS]\nLAMBDA = 0.6 W/mK\nMEW = 15 -\n[MOISTURE_STORAGE]\nFUNCTION = Theta_l(pC)_de\n1 5 9\n0.34 0.17 0.001\n[MOISTURE_TRANSPORT]\nFUNCTION = lgKl(Theta_l)\n0.001 0.17 0.34\n-20 -14 -9\n')
        m=read_material(p);f=m.fields(np.array([.17]),np.array([25.]))
        self.assertAlmostEqual(f[2][0],1e-14,delta=1e-27)
        self.assertAlmostEqual(f[0][0],-1e5,delta=1e-5)
    def test_m6_forward_and_inverse_retention_pair_accepted(self):
        p=self.root/'DUAL_RETENTION_FIXTURE.m6'
        p.write_text('[STORAGE_BASE_PARAMETERS]\nRHO = 1800 kg/m3\nCE = 850 J/kgK\nTHETA_POR = 0.38 m3/m3\nTHETA_EFF = 0.34 m3/m3\n[TRANSPORT_BASE_PARAMETERS]\nLAMBDA = 0.6 W/mK\nMEW = 15 -\n[MOISTURE_STORAGE]\nFUNCTION = Theta_l(pC)_de\n1 5 9\n0.34 0.17 0.001\nFUNCTION = pC(Theta_l)_de\n0.001 0.17 0.34\n9 5 1\n[MOISTURE_TRANSPORT]\nFUNCTION = lgKl(Theta_l)\n0.001 0.17 0.34\n-20 -14 -9\n')
        m=read_material(p)
        self.assertEqual(m.store,'Theta_l(pC)_de')
        self.assertAlmostEqual(m.pressure(np.array([.17]),25.)[0],-1e5,delta=1e-5)
    def test_storage_plateau_rejected(self):
        o=copy.deepcopy(self.obj);k=next(k for k in o['functions'] if 'pC' in k);o['functions'][k][1][2]=o['functions'][k][1][1]
        with self.assertRaises(InputError):Material(o)
    def test_wrong_inverse_pressure_direction_rejected(self):
        o=copy.deepcopy(self.obj);o['functions']={k:v for k,v in o['functions'].items() if 'pC' not in k}
        o['functions']['pC(Theta_l)_de']=[[.001,.17,.34],[1,5,9]]
        with self.assertRaises(InputError):Material(o)
    def test_scalars_cannot_replace_liquid_curve(self):
        o=copy.deepcopy(self.obj);o['functions'].pop('lgKl(Theta_l)')
        with self.assertRaises(InputError):Material(o)
    def test_unknown_function_rejected(self):
        o=copy.deepcopy(self.obj);o['functions']['UNSUPPORTED']=[[0,1],[1,2]]
        with self.assertRaises(InputError):Material(o)
    def test_binary_source_rejected(self):
        p=self.root/'binary.m6';p.write_bytes(b'BINARY\x00M6')
        with self.assertRaises(InputError):read_material(p)
    def test_thermal_table_no_extrapolation(self):
        o=copy.deepcopy(self.obj);o['functions']['lambda(T)']=[[10,30],[.5,.6]];m=Material(o)
        with self.assertRaises(InputError):m.fields(np.array([.17]),np.array([40]))

class PhysicsTests(Sandbox):
    def test_closed_uniform_state_is_stationary(self):
        c=read_climate(self.climate());cfg=self.cfg;cfg['_closed_test']=True;cfg['_zero_gravity_test']=True
        m=ColumnModel(cfg,Material(self.obj));state=m.initial();state.sulfate[:]=.1
        final,h,_,_=m.integrate(c,state)
        np.testing.assert_allclose(final.matrix(),state.matrix(),rtol=1e-10,atol=1e-10)
        self.assertEqual(len(h),24);self.assertEqual(h.elapsed_hours.iloc[-1],24)
    def test_water_and_sulfate_mass_closure(self):
        p=Project(ROOT/'projects/demo');c=p.climate('DEMO').subset(0,72);m=ColumnModel(self.cfg,Material(self.obj))
        _,h,_,stats=m.integrate(c)
        self.assertLess(stats['water_max_closure_error_kg_m2'],1e-8)
        self.assertLess(stats['sulfate_max_closure_error_mol_m2'],1e-10)
        self.assertGreater(h.sulfate_total_mol_m2_footprint.iloc[-1],0)
    def test_fixed_support_sulfate_concentration_metric_reported(self):
        p=Project(ROOT/'projects/demo');c=p.climate('DEMO').subset(0,24);m=ColumnModel(self.cfg,Material(self.obj))
        _,h,_,_=m.integrate(c)
        self.assertIn('sulfate_cmax_0p10m_bandavg_mol_m3_liquid',h.columns)
        self.assertTrue(np.isfinite(h.sulfate_cmax_0p10m_bandavg_mol_m3_liquid).all())
        self.assertTrue((h.sulfate_cmax_0p10m_bandavg_mol_m3_liquid <= h.sulfate_cmax_mol_m3_liquid + 1e-12).all())
    def test_discrete_diffusion_eigenmode(self):
        cfg=self.cfg;cfg['_closed_test']=True;cfg['_zero_gravity_test']=True;cfg['salt']['D_sulfate_m2_s']=1e-6
        material=Material(self.obj)
        def constant_fields(theta,T):
            a=np.ones_like(theta);return a*0,a*.5,a*1e-30,a*1e-30,a*.6
        material.fields=constant_fields
        model=ColumnModel(cfg,material);theta=np.full(model.n,.17)
        n0=.17*(2+.5*np.cos(np.pi*model.z/model.H))
        state=State(theta,np.full(model.n,25.),n0)
        final,_,_,_=model.integrate(read_climate(self.climate()),state)
        rate=4e-6*np.sin(np.pi/(2*model.n))**2/model.dz**2
        exact=.17*(2+.5*np.cos(np.pi*model.z/model.H)*np.exp(-rate*24*3600))
        np.testing.assert_allclose(final.sulfate,exact,rtol=3e-6,atol=1e-7)
        self.assertAlmostEqual(final.sulfate.sum(),n0.sum(),places=10)
    def test_solubility_partition_conserves_total(self):
        p=self.root/'phase.csv';pd.DataFrame({'T_C':[0,60],'Na2SO4_mol_kg_water':[.01,.02],'source':['SYNTHETIC TEST ONLY']*2}).to_csv(p,index=False)
        phase=Solubility(p);m=ColumnModel(self.cfg,Material(self.obj),phase)
        th=np.array([.2,.2]);T=np.array([25.,25.]);nt=np.array([10.,.1]);c,solid=m.partition(th,T,nt)
        np.testing.assert_allclose(th*c+solid,nt);self.assertGreater(solid[0],0)
        with self.assertRaises(InputError):phase.concentration(np.array([70.]))
    def test_solubility_without_source_rejected(self):
        p=self.root/'phase.csv';pd.DataFrame({'T_C':[0,60],'Na2SO4_mol_kg_water':[.01,.02],'source':['','']}).to_csv(p,index=False)
        with self.assertRaises(InputError):Solubility(p)
    def test_subzero_forcing_blocked(self):
        m=ColumnModel(self.cfg,Material(self.obj))
        with self.assertRaises(InputError):m.integrate(read_climate(self.climate(tas_C=-1)))

class WorkflowTests(Sandbox):
    def test_real_site_does_not_accept_synthetic_material(self):
        p=self.demo('real');c=p.config;c['is_demo']=False;atomic_json(p.path/'project.json',c);p=Project(p.path)
        with self.assertRaises(InputError):p.material()
    def test_missing_material_blocks_Mohenjo_transport(self):
        p=Project(ROOT/'projects/mohenjo_daro')
        p.config['material_file']='materials/DELIBERATELY_ABSENT_FOR_TEST.m6'
        self.assertFalse(p.readiness()['trial_ready'])
    def test_settings_hash_changes(self):
        p=self.demo('hash');before=p.physics_hash();p.config['boundary']['capillary_supply_mm_day']+=.1
        self.assertNotEqual(before,p.physics_hash())
    def test_invalid_numerical_fields_report_errors(self):
        p=self.demo('settings');p.config['numerics']['spinup_max_cycles']=0;p.config['geometry']['cells']=None
        self.assertTrue(p.validate_settings())
    def test_freeze_blocks_without_qa_and_spinup(self):
        with self.assertRaises(InputError):freeze(self.demo('freeze'))
    def test_cancel_resume_matches_uninterrupted(self):
        interrupted=self.demo('interrupted');uninterrupted=self.demo('uninterrupted');flag=[False]
        def progress(f,msg):
            if f>=1/3:flag[0]=True
        with self.assertRaises(Cancelled):run_transport(interrupted,'DEMO',max_hours=72,progress=progress,cancel=lambda:flag[0])
        cps=list((interrupted.path/'outputs').rglob('checkpoint.json'));self.assertEqual(len(cps),1)
        checkpoint=json.loads(cps[0].read_text());self.assertEqual(checkpoint['completed_hours'],24)
        self.assertTrue((cps[0].parent/checkpoint['state_file']).exists())
        a=run_transport(interrupted,'DEMO',max_hours=72);b=run_transport(uninterrupted,'DEMO',max_hours=72)
        da=pd.read_csv(interrupted.path/a/'hourly_results.csv.gz');db=pd.read_csv(uninterrupted.path/b/'hourly_results.csv.gz')
        self.assertEqual(len(da),72)
        np.testing.assert_allclose(da[['theta_mean_m3_m3','T_mean_C','sulfate_total_mol_m2_footprint']],db[['theta_mean_m3_m3','T_mean_C','sulfate_total_mol_m2_footprint']],rtol=1e-7,atol=1e-8)
        self.assertEqual(run_transport(interrupted,'DEMO',max_hours=72),a)
    def test_closed_uniform_QA_and_spinup_pass_but_demo_freeze_refused(self):
        p=self.demo('checks');c=p.config;c['_closed_test']=True;c['_zero_gravity_test']=True;p.save(c)
        qa=numerical_qa(p,'DEMO',hours=24);self.assertTrue(qa['passed'])
        result=spinup(p,'DEMO');self.assertTrue(result['passed']);self.assertTrue((p.path/'spinup.npz').exists())
        self.assertTrue(p.readiness()['qa_valid']);self.assertTrue(p.readiness()['spinup_valid'])
        with self.assertRaisesRegex(InputError,'Synthetic'):freeze(p)
    def test_QA_record_invalid_after_climate_change(self):
        p=self.demo('qa_stale');c=p.config;c['_closed_test']=True;c['_zero_gravity_test']=True;p.save(c)
        numerical_qa(p,'DEMO',hours=24);self.assertTrue(p.readiness()['qa_valid'])
        file=p.file(p.run('DEMO')['file']);file.write_bytes(file.read_bytes()+b'\n')
        self.assertFalse(p.readiness()['qa_valid'])
    def test_new_site_does_not_inherit_salt_or_supply(self):
        shutil.copytree(ROOT/'projects/mohenjo_daro',self.root/'mohenjo_daro',ignore=shutil.ignore_patterns('climate','outputs'))
        sid=new_project(self.root,'New site');p=Project(self.root/sid)
        self.assertIsNone(p.config['salt']['D_sulfate_m2_s']);self.assertIsNone(p.config['boundary']['capillary_supply_mm_day'])
        self.assertFalse(p.config['runs']);self.assertFalse(p.readiness()['trial_ready'])

if __name__=='__main__':unittest.main(verbosity=2)
