"""Build any comparison table straight from the eval cache (no text scraping).

eval_round10_A.py writes scripts/analysis/.eval_cache/<run>.json holding the 16
scalars per run.  Reading those is exact and instant, and avoids the fixed-width
parsing that silently mis-aligned columns when new runs with different label
widths were added.

Usage:  python3 table.py [run_key ...]      default: a standard comparison set
"""
import json, os, sys

CD = os.path.join(os.path.dirname(os.path.abspath(__file__)), '.eval_cache')
ROWS = [('SOUTHERN OCEAN 45-65S ocean', None),
        ('  TOA SW CRE [W/m2]', 'so_cre'), ('  cloud area [%]', 'so_cld'),
        ('SIBERIA land JJA', None),
        ('  2m temperature [C]', 'sib_t'), ('  surface net SW', 'sib_sw'),
        ('  cloud area [%]', 'sib_cld'),
        ('GLOBAL guardrails', None),
        ('  net TOA [W/m2]', 'g_toa'), ('  surface flux [W/m2]', 'g_sfc'),
        ('  tropics net TOA', 'trop'),
        ('SW RMSE vs CERES (ocean)', None),
        ('  Southern Ocean', 'rmse_so'), ('  subpolar N Atl', 'rmse_spna'),
        ('  Nordic Seas', 'rmse_nordic'), ('  global SW', 'rmse_sw'),
        ('  T2m vs ERA5 [K]', 'rmse_t2m')]
DEFAULT = ['amip_pi_base', 'amip_picontrol', 'amip_B5_capdcycl0', 'amip_D2b_inpsea02_p700',
           'amip_F1_z0h10', 'amip_F2_lai3', 'amip_F3_cov07', 'amip_F4_rsmin1000',
           'amip_F5_allveg']
runs = sys.argv[1:] or DEFAULT
data = {}
for r in runs:
    p = os.path.join(CD, r + '.json')
    if os.path.exists(p):
        data[r] = json.load(open(p))['m']
    else:
        print(f'  !! no cache for {r} -- run eval_round10_A.py first')
short = {r: r.replace('amip_', '')[:11] for r in data}
w = 12
print('\n' + ' ' * 26 + ''.join(f'{short[r]:>{w}s}' for r in data))
for label, key in ROWS:
    if key is None:
        print(f'-- {label}')
        continue
    print(f'{label:26s}' + ''.join(f'{data[r][key]:{w}.2f}' for r in data))
