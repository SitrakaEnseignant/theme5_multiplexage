"""Exécution réelle du graphe GNU Radio du Thème 5 sans écran, mesures et figures.
Les afficheurs QT GUI sont remplacés par des Vector Sink et le curseur QT GUI Range par une variable ;
aucun bloc de traitement n'est modifié, sauf le bloc Throttle, réglé à une cadence très élevée.
    python3 En_ligne/executer_grc_ci.py"""
import os, sys, json, time, threading, subprocess, importlib, yaml
import numpy as np
NOMF = {'bbag': 'bruit seul', 'echo_court': 'écho court', 'echo_long': 'écho long'}
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
ICI = os.path.dirname(os.path.abspath(__file__)); RACINE = os.path.dirname(ICI)
sys.path.insert(0, ICI)
from ofdm import FS, N, CP, NS, DF, NPRE, NSYM, OCC, ND, QPSK, LAB, CANAUX, sigma2, H, teb_theorique, teb_canal_connu
G = os.path.join(RACINE, 'GNU_Radio'); OUT = os.path.join(RACINE, 'resultats_gnuradio'); HL = os.path.join(OUT, 'hl')
os.makedirs(HL, exist_ok=True)
SINKS = ('qtgui_time_sink_x', 'qtgui_freq_sink_x', 'qtgui_const_sink_x', 'qtgui_eye_sink_x')
EBN0 = {'bbag': [0, 2, 4, 6, 8, 10], 'echo_court': [0, 4, 8, 12, 16, 20]}


def headless(grc, valeurs, params, tag):
    g = yaml.safe_load(open(grc, encoding='utf8'))
    fid = g['options']['parameters']['id'] + '_hl_' + tag
    g['options']['parameters'].update({'id': fid, 'generate_options': 'no_gui', 'run_options': 'run'})
    blocs, sinks = [], {}
    for b in g['blocks']:
        p = b['parameters']
        if b['id'] == 'qtgui_range':
            b = {'name': b['name'], 'id': 'variable', 'parameters': {'comment': '', 'value': p['value']}, 'states': b['states']}
        if b['name'] in valeurs: b['parameters']['value'] = str(valeurs[b['name']])
        if b['name'] in params: b['parameters'].update({k: str(v) for k, v in params[b['name']].items()})
        if b['id'] == 'blocks_throttle': b['parameters']['samples_per_second'] = '1e9'
        if b['id'] in SINKS:
            n = int(p.get('nconnections', 1)); typ = p['type']; sinks[b['name']] = (typ, n)
            for k in range(n):
                blocs.append({'name': f"{b['name']}_{k}", 'id': 'blocks_vector_sink_x', 'states': b['states'],
                              'parameters': {'affinity': '', 'alias': '', 'comment': '', 'type': typ, 'vlen': '1',
                                             'reserve_items': '1024'}})
            continue
        blocs.append(b)
    conns = [[c[0], c[1], f"{c[2]}_{c[3]}", '0'] if c[2] in sinks else c for c in g['connections']]
    g['blocks'], g['connections'] = blocs, conns
    path = os.path.join(HL, fid + '.grc'); yaml.safe_dump(g, open(path, 'w', encoding='utf8'), allow_unicode=True, sort_keys=False)
    return path, fid, sinks


def executer(ntrames, taps, ebn0=None, tag='a'):
    amp = 0.0 if ebn0 is None else float(np.sqrt(sigma2(10 ** (ebn0 / 10)) / 2))
    path, fid, sinks = headless(os.path.join(G, 'ofdm_tx_rx_sim.grc'), {'ntrames': ntrames, 'noise_amp': amp, 'taps': list(taps)},
                                {'src': {'repeat': 'False'}}, tag)
    r = subprocess.run(['grcc', '-o', HL, path], capture_output=True, text=True)
    if r.returncode: raise RuntimeError(r.stdout + r.stderr + '\n' + diagnostic(path))
    sys.path.insert(0, HL); mod = importlib.reload(importlib.import_module(fid))
    tb = getattr(mod, fid)(); t0 = time.time(); tb.start()
    th = threading.Thread(target=tb.wait); th.start(); th.join(timeout=300)
    fini = not th.is_alive()
    if not fini: tb.stop(); th.join()
    d = {s: [np.array(getattr(tb, f'{s}_{k}').data()) for k in range(n)] for s, (typ, n) in sinks.items()}
    return d, {'duree_s': round(time.time() - t0, 2), 'fin_naturelle': fini, 'ecart_type_bruit': amp}


def diagnostic(path):
    """en cas d'échec de grcc : nom du bloc concerné par chaque erreur de validation, ports du bloc Python"""
    try:
        from gnuradio import gr
        from gnuradio.grc.core.platform import Platform
        p = Platform(name='diagnostic', prefs=gr.prefs(), version=gr.version(),
                     version_parts=(gr.major_version(), gr.api_version(), gr.minor_version()))
        p.build_library(); fg = p.make_flow_graph(); fg.grc_file_path = path
        fg.import_data(p.parse_flow_graph(path)); fg.rewrite(); fg.validate()
        out = []
        for elem, msg in fg.iter_error_messages():
            b = getattr(elem, 'parent_block', elem); out.append(f'BLOC {getattr(b, "name", "?")} / {elem} : {msg}')
        for b in fg.blocks:
            if b.key == 'epy_block':
                out.append('ports du bloc Python : ' + str([(pp.name, pp.dtype, pp.vlen) for pp in list(b.sinks) + list(b.sources)]))
        return '\n'.join(out)
    except Exception as e:
        return f'diagnostic impossible : {e!r}'


def erreurs(d):
    tx = d['symboles'][0].astype(int); dec = d['symboles'][1].astype(int); n = min(len(tx), len(dec))
    lab = np.array(LAB); eb = int(sum(bin(int(v)).count('1') for v in lab[tx[:n]] ^ lab[dec[:n]]))
    return int(np.sum(tx[:n] != dec[:n])), n, eb, 2 * n


def evm(d):
    tx = d['symboles'][0].astype(int); z = d['constellation_rx'][0]; n = min(len(tx), len(z))
    ref = np.array(QPSK)[tx[:n]]; return float(np.sqrt(np.mean(np.abs(z[:n] - ref) ** 2)))


def dsp(x, Nf=1024):
    w = np.hanning(Nf); P = np.mean([np.abs(np.fft.fftshift(np.fft.fft(x[i:i + Nf] * w))) ** 2 for i in range(0, len(x) - Nf, Nf // 2)], 0)
    return (np.arange(Nf) - Nf // 2) * FS / Nf, P


R, ECHECS = {}, []
# 1. sans bruit, trois canaux
for nom, taps in CANAUX.items():
    d, info = executer(10, taps, None, 'nb_' + nom)
    es, ns, eb, nb = erreurs(d); e = evm(d)
    hest = d['canal_rx'][0][:N] / N; hth = np.array([np.sum(np.asarray(taps) * np.exp(-2j * np.pi * k * np.arange(len(taps)) / N))
                                                   for k in range(-N // 2, N // 2)])
    occ = np.array(OCC) + N // 2; eh = float(np.max(np.abs(hest[occ] - hth[occ])))
    R.setdefault(nom, {})['sans_bruit'] = {'erreurs_symboles': es, 'symboles': ns, 'TEB': eb / nb, 'evm': round(e, 6),
                                           'ecart_max_canal_estime': round(eh, 6), **info}
    print(nom, R[nom]['sans_bruit'], flush=True)
    if nom == 'echo_long':
        if e < 0.1: ECHECS.append(f'echo_long : EVM {e:.3f} alors qu\'une dégradation est attendue')
    else:
        if es > 0 or e > 1e-3 or eh > 1e-3: ECHECS.append(f'{nom} : {es} erreurs, EVM {e:.2e}, écart canal {eh:.2e} sans bruit')
    R[nom]['donnees_sans_bruit'] = d
# spectre, puissance, PAPR sur le signal émis sans bruit (canal idéal)
s = R['bbag']['donnees_sans_bruit']['scope_emis'][0]
f, P = dsp(s); frac = float(P[np.abs(f) <= 26.5 * DF].sum() / P.sum()); pm = float(np.mean(np.abs(s) ** 2))
X = s[:len(s) // NS * NS].reshape(-1, NS)[:, CP:]; papr = 10 * np.log10(np.max(np.abs(X) ** 2, 1) / pm)
p8 = float(np.mean(papr > 8)); p8t = float(1 - (1 - np.exp(-10 ** 0.8)) ** 52)
R['signal'] = {'part_puissance_bande': round(frac, 5), 'puissance_moyenne': round(pm, 3), 'puissance_attendue': 52.0,
               'papr_moyen_dB': round(float(papr.mean()), 2), 'P_papr_sup_8dB': round(p8, 4), 'P_papr_sup_8dB_approx': round(p8t, 4)}
print(R['signal'], flush=True)
if frac < 0.98: ECHECS.append(f'spectre : {frac:.3f} de la puissance dans la bande')
if abs(pm / 52 - 1) > 0.05: ECHECS.append(f'puissance moyenne {pm:.2f} au lieu de 52')
# 2. taux d'erreur binaire
for nom, ebs in EBN0.items():
    R[nom]['taux'] = []
    for e_db in ebs:
        g = 10 ** (e_db / 10); d, info = executer(100, CANAUX[nom], e_db, f'{nom}_{e_db}')
        es, ns, eb, nb = erreurs(d); th = teb_theorique(CANAUX[nom], g); mes = eb / nb
        ok = eb < 50 or 1 / 1.5 < mes / th < 1.5
        R[nom]['taux'].append({'EbN0_dB': e_db, 'TEB': mes, 'erreurs': eb, 'bits': nb, 'theorie': th,
                               'canal_connu': teb_canal_connu(CANAUX[nom], g), 'conforme': bool(ok), **info})
        print(nom, e_db, f'{mes:.3e}', eb, f'{th:.3e}', ok, flush=True)
        if not ok: ECHECS.append(f'{nom} : Eb/N0 = {e_db} dB, TEB {mes:.3g}, théorie {th:.3g}')
        if nom == 'echo_court' and e_db == 20: R[nom]['donnees_20dB'] = d
# 3. figures
fig, ax = plt.subplots(2, 3, figsize=(17, 10))
t = np.arange(3 * NS) / FS * 1e6; s0 = s[2 * NS:5 * NS]
ax[0, 0].plot(t, s0.real, lw=0.8, label='partie réelle'); ax[0, 0].plot(t, np.abs(s0), lw=0.8, label='module')
for k in range(3): ax[0, 0].axvspan(k * NS, k * NS + CP, color='orange', alpha=0.25)
ax[0, 0].set_xlabel('temps (µs)'); ax[0, 0].set_title('signal émis (scope-emis) ; préfixes cycliques en orange'); ax[0, 0].legend(fontsize=8)
ax[0, 1].plot(f / 1e3, 10 * np.log10(P / P.max() + 1e-12)); ax[0, 1].set_ylim(-60, 3); ax[0, 1].set_xlabel('fréquence (kHz)')
ax[0, 1].set_title('spectre du signal émis (spectre-emis)'); ax[0, 1].set_ylabel('dB')
z = ax[0, 2]; zz = np.sort(papr); z.semilogy(zz, 1 - np.arange(len(zz)) / len(zz), '.', ms=3, label='mesure')
zt = np.linspace(4, 12, 50); z.semilogy(zt, 1 - (1 - np.exp(-10 ** (zt / 10))) ** 52, label=r'approximation $1-(1-e^{-z})^{52}$')
z.set_xlabel('PAPR (dB)'); z.set_ylabel('P(PAPR > z)'); z.set_title('rapport puissance crête sur puissance moyenne'); z.legend(fontsize=8); z.set_ylim(1e-3, 1.2)
d20 = R['echo_court']['donnees_20dB']; br = d20['constellation_brute'][0][:N * 400].reshape(-1, N)[:, np.array(OCC) + N // 2].ravel() / N
ax[1, 0].plot(br.real, br.imag, '.', ms=2); ax[1, 0].set_aspect('equal'); ax[1, 0].set_title('porteuses avant égalisation, écho court, 20 dB')
zr = d20['constellation_rx'][0][:20000]; ax[1, 1].plot(zr.real, zr.imag, '.', ms=2); ax[1, 1].set_aspect('equal')
ax[1, 1].set_xlim(-2, 2); ax[1, 1].set_ylim(-2, 2); ax[1, 1].set_title('après égalisation (constellation-rx), écho court, 20 dB')
kk = np.arange(-N // 2, N // 2)
for nom, c in (('echo_court', 'C0'), ('echo_long', 'C3')):
    he = R[nom]['donnees_sans_bruit']['canal_rx'][0][:N] / N; ht = np.abs([np.sum(np.asarray(CANAUX[nom]) * np.exp(-2j * np.pi * k * np.arange(len(CANAUX[nom])) / N)) for k in kk])
    ax[1, 2].plot(kk[np.array(OCC) + 32], np.abs(he[np.array(OCC) + 32]), 'o', ms=3, color=c, label=f'{NOMF[nom]} estimé')
    ax[1, 2].plot(kk, ht, '-', lw=0.8, color=c, label=f'{NOMF[nom]} théorique')
ax[1, 2].set_xlabel('indice de porteuse'); ax[1, 2].set_title('module du canal estimé (canal-rx), sans bruit'); ax[1, 2].legend(fontsize=7)
for a in ax.flat: a.grid(ls=':')
fig.suptitle('ofdm-tx-rx-sim - GNU Radio'); fig.tight_layout(); fig.savefig(os.path.join(OUT, 'ofdm_tx_rx_sim.png'), dpi=110); plt.close(fig)
fig, a = plt.subplots(figsize=(9, 6))
for i, nom in enumerate(EBN0):
    T = R[nom]['taux']; e = np.array([p['EbN0_dB'] for p in T]); y = np.array([p['TEB'] for p in T]); ee = np.linspace(e.min(), e.max(), 40)
    a.semilogy(ee, [teb_theorique(CANAUX[nom], 10 ** (x / 10)) for x in ee], color=f'C{i}', lw=1, label=f'{NOMF[nom]} : théorie, canal estimé')
    a.semilogy(ee, [teb_canal_connu(CANAUX[nom], 10 ** (x / 10)) for x in ee], color=f'C{i}', lw=1, ls='--', label=f'{NOMF[nom]} : canal connu')
    a.semilogy(e, np.where(y > 0, y, np.nan), 'o', color=f'C{i}', label=f'{NOMF[nom]} : mesure GNU Radio')
a.set_xlabel('Eb/N0 (dB)'); a.set_ylabel('TEB'); a.set_ylim(1e-5, 0.5); a.grid(which='both', ls=':'); a.legend(fontsize=8)
a.set_title('TEB de la QPSK OFDM mesuré dans GNU Radio et théorie'); fig.tight_layout(); fig.savefig(os.path.join(OUT, 'teb_gnuradio.png'), dpi=110); plt.close(fig)
for nom in R:
    R[nom].pop('donnees_sans_bruit', None); R[nom].pop('donnees_20dB', None)
R['echecs'] = ECHECS
json.dump(R, open(os.path.join(OUT, 'mesures.json'), 'w'), indent=1, ensure_ascii=False, default=float)
print('ECHECS :', ECHECS)
sys.exit(1 if ECHECS else 0)
