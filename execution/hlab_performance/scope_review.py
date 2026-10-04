"""Explicit acknowledgement of the equilibrium water-phase extension."""
from pathlib import Path
from .safety import atomic, now, hash_file, read_json, safe_path, SafetyError
POLICY='equilibrium-freeze-thaw-1.0'
TEXT=('Cold cases may continue with a local-equilibrium water/ice partition derived from the frozen material retention curve and the pure-water ice/liquid equilibrium relation. The extension adds fusion latent heat, ice sensible heat, liquid-only transport and salt exclusion from ice while preserving the original arithmetic whenever no ice is present. It is not a frost-damage, snow, nucleation-hysteresis, pore-expansion, mixed-electrolyte, or salt-freezing-point model. The pure-water assumption is conservative with respect to ice formation because salt freezing-point depression is not credited. The phase formulation is declared for material temperatures down to -20 C. Existing completed outputs are retained and not recalculated; this approval is separate from the original numerical QA.')

def policy_files():
    base=Path(__file__).parent
    return {n:hash_file(base/n) for n in ('phase_change.py','pore_ice.py','worker.py')}

def review_path(root,site):
    return safe_path(Path(root)/'performance_runtime'/'scope_review',site+'.json')

def reviewed(root,site,physics):
    try:
        data=read_json(review_path(root,site))
        return (data.get('approved') is True and data.get('policy')==POLICY and data.get('physics_hash')==physics and data.get('policy_files')==policy_files())
    except (OSError,ValueError,KeyError):return False

def require_review(root,site,physics,approve=False):
    if reviewed(root,site,physics):return read_json(review_path(root,site))
    if approve is not True:
        raise SafetyError('Review and tick the equilibrium freeze/thaw extension acknowledgement before running. Existing completed outputs and checkpoints were not changed.')
    receipt={'approved':True,'policy':POLICY,'physics_hash':physics,'policy_files':policy_files(),'approved_utc':now(),'text':TEXT}
    atomic(review_path(root,site),receipt);return receipt
