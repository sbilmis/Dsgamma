"""Prepare (do not apply) the consistent Python/Wolfram source correction."""
import difflib
import json
from pathlib import Path
import sys
import types

import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0,str(HERE))
from candidate_dimension5 import Inputs, heavy_line_mixed_component
from twopoint_ds1_matrix_sumrule import local_components

changes=[]
python_path=ROOT/"scripts/twopoint_ds1_matrix_sumrule.py"
old=python_path.read_text()
marker='    return {\n        "d3_ms0": exp_factor * q0,'
addition='''    # Second d=5 topology: one background gluon on the heavy line,
    # contracted with the local open quark-gluon operator. Its AA entry
    # vanishes; include it in d5_mixed so OPE diagnostics use the total.
    # Derivation: review/mixed_condensate_resolution/.
    d5 += mixed * (1.0 - mc * mc / m2) * np.array(
        [[0.0, 1.0 / (12.0 * den)],
         [1.0 / (12.0 * den), mc / (6.0 * den**2)]], dtype=float,
    )

'''
assert old.count(marker)==1
new=old.replace(marker,addition+marker)
new=new.replace('  * the local mixed condensate <sbar g sigma.G s> (dimension five).',
    '  * both tree-level topologies of the mixed condensate (dimension five).')
new=new.replace('The AA entry reproduces the Colangelo two-point sum rule analytically.  The\n'
    'AB and BB local terms are generated from the same coordinate-space propagator\n'
    'and were independently traced in ``twopoint_local_condensate_audit.wl``.',
    'The AA entry reproduces the Colangelo two-point sum rule analytically. The\n'
    'local terms include the light-bilocal expansion and the heavy-line gluon\n'
    'contribution to the mixed operator. Their complete leading-d5 derivation\n'
    'and independent checks are in ``review/mixed_condensate_resolution/``.')
new=new.replace('    """Borel-transformed d=3 and d=5 local contributions.\n',
    '    """Borel-transformed d=3 and complete leading d=5 local contributions.\n')
new=new.replace('    # Dimension-five mixed condensate from the x^2 propagator term.',
    '    # Light-bilocal part of d=5 from the x^2 propagator term.')
changes.append((python_path,old,new))

for relative,mass,indent in [
    ("scripts/mathematica_twopoint_ds1_matrix_check.wl","mc","  "),
    ("scripts/mathematica_twopoint_bs1_matrix_check.wl","mb","  "),
    ("scripts/build_twopoint_ds1_matrix_notebook.wl","mc","        "),
]:
    path=ROOT/relative;old_wl=path.read_text()
    marker_wl=indent+'<|"d3_ms0" -> e q0, "d3_ms1" -> e q1,'
    assert old_wl.count(marker_wl)==1
    addition_wl=(indent+'(* Heavy-line background-gluon contribution to the same C5 operator. *)\n'
        +indent+f'd5 += mixedSS (1 - {mass}^2/M2) '+ '{{0, 1/(12 den)},\n'
        +indent+f'  {{1/(12 den), {mass}/(6 den^2)}}}};\n')
    new_wl=old_wl.replace(marker_wl,addition_wl+marker_wl)
    if 'notebook' in relative:
        new_wl=new_wl.replace(
            '"- x^2 <sbar g sigma.G s>/192. They are displayed here through their " <>\n'
            '      "Borel-transformed matrix contributions, which are what enter the sum rule."',
            '"- x^2 <sbar g sigma.G s>/192. The leading mixed-condensate matrix " <>\n'
            '      "also includes the heavy-line background-gluon contribution. Both are " <>\n'
            '      "included in d5_mixed below; see review/mixed_condensate_resolution/."')
    changes.append((path,old_wl,new_wl))

patch=''.join(''.join(difflib.unified_diff(old.splitlines(True),new.splitlines(True),
            fromfile='a/'+str(path.relative_to(ROOT)),tofile='b/'+str(path.relative_to(ROOT))))
            for path,old,new in changes)
(HERE/"proposed_implementation.patch").write_text(patch)
(HERE/"proposed_wolfram_sources.json").write_text(json.dumps(
    {str(path.relative_to(ROOT)): new for path,_,new in changes[1:]},indent=2)+"\n")

# Verify the proposed Python edit in memory; never apply it to production.
module=types.ModuleType('d5_patch_check')
module.__file__=str(python_path)
sys.modules[module.__name__]=module
exec(compile(changes[0][2],str(python_path),'exec'),module.__dict__)
errors=[]
for m,ms,K in [(1.27,0,3),(1.27,.126,4.5),(4.18,0,9),(4.18,.126,11)]:
    inp=Inputs(mc=m,ms=ms,qq=-.0165,kappa_s=.8,m0_sq=.8)
    old_parts=local_components(K,inp)
    patched=module.local_components(K,inp)
    expected=old_parts['d5_mixed']+heavy_line_mixed_component(K,inp)
    error=float(np.max(np.abs(patched['d5_mixed']-expected)))
    assert error<1e-15
    for key in ['d3_ms0','d3_ms1','d3_ms2']:
        assert np.array_equal(old_parts[key],patched[key])
    errors.append(error)
(HERE/"proposed_patch_check.json").write_text(json.dumps(dict(
    applied_to_production=False,patched_python_points=len(errors),
    maximum_absolute_error=max(errors),other_local_components_unchanged=True,
    files=[str(p.relative_to(ROOT)) for p,_,_ in changes]),indent=2)+"\n")
print('Prepared patch for',len(changes),'source files; verified Python edit at',len(errors),'points.')
