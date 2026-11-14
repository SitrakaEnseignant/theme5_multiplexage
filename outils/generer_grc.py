"""Génère le graphe de flux GNU Radio Companion 3.10 du Thème 5 (format YAML), en simulation pure :
chaîne OFDM complète (allocation des porteuses, IFFT, préfixe cyclique, canal à trajets multiples et bruit,
retrait du préfixe, FFT, estimation du canal et égalisation, décision).
    python3 outils/generer_grc.py"""
import os, sys, json, copy, yaml
ICI = os.path.dirname(os.path.abspath(__file__)); RACINE = os.path.dirname(ICI)
sys.path.insert(0, os.path.join(RACINE, 'En_ligne'))
from ofdm import FS, N, CP, NPRE, NSYM, OCC, PIL, PILS, ND, QPSK, preambule_decale, table
MODELES = json.load(open(os.path.join(ICI, 'modeles_blocs.json'), encoding='utf8'))
OUT = os.path.join(RACINE, 'GNU_Radio')
COMMUNS = {'alias': '', 'comment': '', 'affinity': '', 'maxoutbuf': '0', 'minoutbuf': '0'}


class Q(str): pass


class FG:
    def __init__(s, fid, titre, desc):
        s.fid, s.titre, s.desc, s.blocks, s.conns = fid, titre, desc, [], []

    def blk(s, nom, bid, x, y, **p):
        base = copy.deepcopy(MODELES.get(bid, COMMUNS))
        for k in [k for k in base if k.startswith(('label', 'color', 'alpha', 'width', 'style', 'marker'))
                  and k[-1].isdigit() and int(''.join(c for c in k if c.isdigit()) or 0) > 10]:
            base.pop(k)
        base.update({k: str(v) for k, v in p.items()})
        base.setdefault('comment', '')
        s.blocks.append({'name': nom, 'id': bid, 'parameters': {k: Q(v) for k, v in base.items()},
                         'states': {'bus_sink': False, 'bus_source': False, 'bus_structure': None,
                                    'coordinate': [int(x), int(y)], 'rotation': 0, 'state': 'enabled'}})
        return nom

    def c(s, a, b, pa=0, pb=0):
        s.conns.append([a, str(pa), b, str(pb)])

    def save(s):
        d = {'options': {'parameters': {'author': 'Thème 5 - multiplexage OFDM (simulation pure)', 'catch_exceptions': 'True',
                 'category': '[GRC Hier Blocks]', 'cmake_opt': '', 'comment': '', 'copyright': '', 'description': s.desc,
                 'gen_cmake': 'On', 'gen_linking': 'dynamic', 'generate_options': 'qt_gui', 'hier_block_src_path': '.:',
                 'id': s.fid, 'max_nouts': '0', 'output_language': 'python', 'placement': '(0,0)', 'qt_qss_theme': '',
                 'realtime_scheduling': '', 'run': 'True', 'run_command': '{python} -u {filename}', 'run_options': 'prompt',
                 'sizing_mode': 'fixed', 'thread_safe_setters': '', 'title': s.titre, 'window_size': '(1300,1000)'},
                 'states': {'bus_sink': False, 'bus_source': False, 'bus_structure': None, 'coordinate': [8, 8],
                            'rotation': 0, 'state': 'enabled'}},
             'blocks': s.blocks, 'connections': s.conns, 'metadata': {'file_format': 1, 'grc_version': '3.10.9.2'}}

        class D(yaml.SafeDumper): pass
        D.add_representer(Q, lambda dp, v: dp.represent_scalar('tag:yaml.org,2002:str', str(v), style="'"))
        with open(os.path.join(OUT, s.fid + '.grc'), 'w', encoding='utf8') as f:
            yaml.dump(d, f, Dumper=D, sort_keys=False, allow_unicode=True, default_flow_style=False, width=1000)


def rng(g, nom, val, a, b, pas, lab, hint, x, y=110):
    return g.blk(nom, 'qtgui_range', x, y, value=val, start=a, stop=b, step=pas, label=lab, gui_hint=hint,
                 widget='counter_slider', rangeType='float', orient='QtCore.Qt.Horizontal', min_len='200', outputmsgname='value')


def time_sink(g, nom, typ, n, labels, hint, x, y, size='1024', ymin='-2', ymax='2'):
    p = dict(type=typ, name=f'"{nom}"', size=size, srate='samp_rate', nconnections=n, ymin=ymin, ymax=ymax, gui_hint=hint,
             autoscale='False', grid='True', legend='True', update_time='0.10', entags='False', ctrlpanel='False', stemplot='False')
    for i, l in enumerate(labels, 1): p[f'label{i}'] = f'"{l}"'
    return g.blk(nom, 'qtgui_time_sink_x', x, y, **p)


def const_sink(g, nom, hint, x, y, lim='2'):
    return g.blk(nom, 'qtgui_const_sink_x', x, y, type='complex', name=f'"{nom}"', size='2048', nconnections='1',
                 xmin=f'-{lim}', xmax=lim, ymin=f'-{lim}', ymax=lim, gui_hint=hint, autoscale='False', grid='True', label1='"porteuses"')


def construire():
    g = FG('ofdm_tx_rx_sim', 'OFDM : émetteur, canal à trajets multiples et récepteur (simulation)',
           "Chaîne OFDM 802.11a simplifiée : 64 porteuses, préfixe de 16 échantillons, 48 porteuses de données en QPSK, "
           "4 pilotes, 2 symboles d'apprentissage par trame. Aucun matériel de radio logicielle.")
    g.blk('imports', 'import', 1500, 12, imports='import numpy as np')
    g.blk('samp_rate', 'variable', 16, 12, value=repr(FS))
    g.blk('fft_len', 'variable', 140, 12, value=str(N))
    g.blk('cp_len', 'variable', 250, 12, value=str(CP))
    g.blk('npre', 'variable', 360, 12, value=str(NPRE))
    g.blk('nsym_trame', 'variable', 470, 12, value=str(NSYM))
    g.blk('ntrames', 'variable', 600, 12, value='100')
    g.blk('occupied_carriers', 'variable', 720, 12, value=str((tuple(OCC),)))
    g.blk('pilot_carriers', 'variable', 900, 12, value=str((tuple(PIL),)))
    g.blk('pilot_symbols', 'variable', 1050, 12, value=str((tuple(PILS),)))
    g.blk('preambule', 'variable', 1200, 12, value=str([float(v.real) for v in preambule_decale()]))
    g.blk('taps', 'variable', 1350, 12, value='[1.0]')
    g.blk('constellation', 'variable_constellation', 1650, 12, const_points=table(QPSK), sym_map='[0, 1, 2, 3]',
          dims='1', normalization='digital.constellation.AMPLITUDE_NORMALIZATION', npwr='1.0', precision='8',
          rot_sym='4', soft_dec_lut='None', type='calcdist')
    rng(g, 'noise_amp', '0.0', '0', '5', '0.01', 'Écart type du bruit par voie (I et Q)', '0,0,1,2', 16)
    y1, y2 = 300, 620
    # émetteur
    g.blk('src', 'analog_random_source_x', 16, y1, type='byte', min='0', max='4', num_samps='ntrames*nsym_trame*48', repeat='True')
    g.blk('symboles_tx', 'digital_chunks_to_symbols_xx', 200, y1, in_type='byte', out_type='complex', symbol_table=table(QPSK),
          dimension='1', num_ports='1', vlen='1')
    g.blk('trame', 'blocks_stream_to_tagged_stream', 380, y1, type='complex', vlen='1', packet_len='nsym_trame*48',
          len_tag_key='"packet_len"')
    g.blk('allocateur', 'digital_ofdm_carrier_allocator_cvc', 560, y1, fft_len='fft_len', occupied_carriers='occupied_carriers',
          pilot_carriers='pilot_carriers', pilot_symbols='pilot_symbols', sync_words='(preambule, preambule)',
          len_tag_key='"packet_len"', output_is_shifted='True')
    g.blk('ifft_tx', 'fft_vxx', 800, y1, type='complex', fft_size='fft_len', forward='False', window='[1.0] * fft_len',
          shift='True', nthreads='1')
    g.blk('prefixe', 'digital_ofdm_cyclic_prefixer', 1000, y1, input_size='fft_len', cp_len='cp_len', rolloff='0',
          tagname='"packet_len"')
    g.blk('throttle', 'blocks_throttle', 1200, y1, type='complex', samples_per_second='samp_rate', ignoretag='True', vlen='1')
    g.c('src', 'symboles_tx'); g.c('symboles_tx', 'trame'); g.c('trame', 'allocateur'); g.c('allocateur', 'ifft_tx')
    g.c('ifft_tx', 'prefixe'); g.c('prefixe', 'throttle')
    # canal : trajets multiples puis bruit gaussien complexe
    g.blk('canal_trajets', 'filter_fir_filter_xxx', 1400, y1, type='ccc', decim='1', taps='taps', samp_delay='0')
    g.blk('bruit_I', 'analog_noise_source_x', 1400, y1 + 110, type='float', noise_type='analog.GR_GAUSSIAN', amp='noise_amp', seed='1')
    g.blk('bruit_Q', 'analog_noise_source_x', 1400, y1 + 200, type='float', noise_type='analog.GR_GAUSSIAN', amp='noise_amp', seed='2')
    g.blk('bruit', 'blocks_float_to_complex', 1600, y1 + 150, vlen='1')
    g.blk('canal', 'blocks_add_xx', 1800, y1, type='complex', num_inputs='2', vlen='1')
    g.c('throttle', 'canal_trajets'); g.c('canal_trajets', 'canal'); g.c('bruit_I', 'bruit'); g.c('bruit_Q', 'bruit', 0, 1)
    g.c('bruit', 'canal', 0, 1)
    time_sink(g, 'scope_emis', 'complex', '1', ['signal émis'], '1,0,1,1', 1200, y1 - 160, size='400', ymin='-30', ymax='30')
    g.blk('spectre_emis', 'qtgui_freq_sink_x', 1400, y1 - 160, type='complex', name='"spectre_emis"', fftsize='1024',
          bw='samp_rate', fc='0', nconnections='1', average='0.05', ymin='-80', ymax='40', gui_hint='1,1,1,1', autoscale='False',
          grid='True', legend='True', freqhalf='True', label1='"signal émis"', wintype='window.WIN_HANN')
    g.c('throttle', 'scope_emis'); g.c('throttle', 'spectre_emis')
    # récepteur
    g.blk('retrait_prefixe', 'blocks_keep_m_in_n', 16, y2, type='complex', m='fft_len', n='fft_len + cp_len', offset='cp_len')
    g.blk('serie_parallele', 'blocks_stream_to_vector', 220, y2, type='complex', num_items='fft_len', vlen='1')
    g.blk('fft_rx', 'fft_vxx', 420, y2, type='complex', fft_size='fft_len', forward='True', window='[1.0] * fft_len',
          shift='True', nthreads='1')
    src = open(os.path.join(ICI, 'egaliseur_epy.py'), encoding='utf8').read()
    g.blk('egaliseur', 'epy_block', 640, y2, _source_code=src, fft_len='fft_len', npre='npre', nsym='nsym_trame',
          preambule='preambule', porteuses=str(list(OCC)))
    g.blk('donnees', 'blocks_vector_to_stream', 880, y2, type='complex', num_items=str(ND), vlen='1')
    g.blk('canal_estime', 'blocks_vector_to_stream', 880, y2 + 120, type='complex', num_items='fft_len', vlen='1')
    g.blk('porteuses_brutes', 'blocks_vector_to_stream', 640, y2 + 200, type='complex', num_items='fft_len', vlen='1')
    g.blk('decision', 'digital_constellation_decoder_cb', 1080, y2, constellation='constellation')
    g.blk('dec_float', 'blocks_char_to_float', 1260, y2, scale='1', vlen='1')
    g.blk('tx_float', 'blocks_char_to_float', 200, y1 - 110, scale='1', vlen='1')
    g.c('canal', 'retrait_prefixe'); g.c('retrait_prefixe', 'serie_parallele'); g.c('serie_parallele', 'fft_rx')
    g.c('fft_rx', 'egaliseur'); g.c('egaliseur', 'donnees', 0, 0); g.c('egaliseur', 'canal_estime', 1, 0)
    g.c('fft_rx', 'porteuses_brutes'); g.c('donnees', 'decision'); g.c('decision', 'dec_float'); g.c('src', 'tx_float')
    const_sink(g, 'constellation_rx', '2,1,1,1', 1080, y2 + 120)
    const_sink(g, 'constellation_brute', '2,0,1,1', 880, y2 + 220, lim='30')
    time_sink(g, 'canal_rx', 'complex', '1', ['canal estimé'], '3,0,1,1', 1080, y2 + 220, size='64', ymin='-3', ymax='3')
    time_sink(g, 'symboles', 'float', '2', ['symbole émis', 'symbole décidé'], '3,1,1,1', 1460, y2, size='96', ymin='-0.5', ymax='3.5')
    g.c('donnees', 'constellation_rx'); g.c('porteuses_brutes', 'constellation_brute'); g.c('canal_estime', 'canal_rx')
    g.c('tx_float', 'symboles', 0, 0); g.c('dec_float', 'symboles', 0, 1)
    g.save()
    print('écrit', g.fid, len(g.blocks), 'blocs', len(g.conns), 'liaisons')


def parametres_json():
    """mêmes paramètres pour MATLAB (lus par jsondecode dans le script Simulink)"""
    from ofdm import FS, N, CP, NPRE, NSYM, OCC, PIL, PILS, USED, QPSK, LAB, CANAUX
    d = {'fs': FS, 'N': N, 'CP': CP, 'NPRE': NPRE, 'NSYM': NSYM, 'OCC': OCC, 'PIL': PIL, 'PILS': PILS, 'USED': USED,
         'PRE': [float(v.real) for v in preambule_decale()], 'QPSK_re': [p.real for p in QPSK], 'QPSK_im': [p.imag for p in QPSK],
         'LAB': LAB, 'canaux': {k: [float(x) for x in v] for k, v in CANAUX.items()}}
    json.dump(d, open(os.path.join(RACINE, 'Simulink', 'parametres_ofdm.json'), 'w'), indent=1)


if __name__ == '__main__':
    construire()
    parametres_json()
