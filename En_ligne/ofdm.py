"""Paramètres communs de la chaîne OFDM du Thème 5 (Bloc 2), repris à l'identique par le script MATLAB
Simulink/construire_et_verifier_modeles.m et par le graphe GNU Radio GNU_Radio/ofdm_tx_rx_sim.grc.
Disposition des porteuses du 802.11a, comme le graphe lesson30_OFDM du cours Learn SDR."""
import numpy as np
from scipy.special import erfc

FS = 1e6                 # cadence d'échantillonnage (Hz), comme la leçon 30 de Learn SDR
N = 64                   # taille de la FFT
CP = 16                  # préfixe cyclique (échantillons), N/4
NS = N + CP              # échantillons par symbole OFDM
DF = FS / N              # écart entre sous-porteuses : 15 625 Hz
NPRE = 2                 # symboles d'apprentissage (préambule) au début de chaque trame
NSYM = 20                # symboles OFDM de données par trame
# porteuses de données et pilotes du 802.11a (indices de -32 à 31, la porteuse 0 est vide)
OCC = list(range(-26, -21)) + list(range(-20, -7)) + list(range(-6, 0)) + list(range(1, 7)) + list(range(8, 21)) + list(range(22, 27))
PIL = [-21, -7, 7, 21]
PILS = [1, 1, 1, -1]
ND = len(OCC)            # 48 porteuses de données
USED = sorted(OCC + PIL) # 52 porteuses utilisées
# préambule : ±1 sur les 52 porteuses utilisées (suite pseudo-aléatoire fixe), 0 ailleurs
_r = np.random.default_rng(802)
PRE = {c: float(v) for c, v in zip(USED, 2 * _r.integers(0, 2, len(USED)) - 1)}
# QPSK, code de Gray : indice m -> point ; voisins à un bit près
QPSK = [complex(1, 1) / np.sqrt(2), complex(-1, 1) / np.sqrt(2), complex(-1, -1) / np.sqrt(2), complex(1, -1) / np.sqrt(2)]
LAB = [0, 1, 3, 2]
M = 4; BITS = 2
# canaux étudiés : réponse impulsionnelle en échantillons (1 µs par échantillon)
CANAUX = {'bbag': [1.0], 'echo_court': [1.0, 0, 0, 0.9], 'echo_long': [1.0] + [0] * 19 + [0.9]}


def preambule_decale():
    """symbole d'apprentissage au format décalé de GNU Radio (porteuse 0 à l'indice 32)"""
    v = np.zeros(N, complex)
    for c, x in PRE.items(): v[c + N // 2] = x
    return v


def sigma2(ebn0_lin):
    """variance totale (I + Q) du bruit complexe par échantillon donnant le Eb/N0 utile voulu.
    Avec l'IFFT et la FFT non normalisées de GNU Radio, le rapport par porteuse vaut N |X|^2 / sigma^2."""
    return N / (BITS * ebn0_lin)


def Q(x): return 0.5 * erfc(np.asarray(x) / np.sqrt(2))


def H(taps):
    """réponse en fréquence du canal sur les porteuses de données"""
    h = np.asarray(taps, complex); k = np.asarray(OCC)
    return np.array([np.sum(h * np.exp(-2j * np.pi * kk * np.arange(len(h)) / N)) for kk in k])


def teb_theorique(taps, ebn0_lin, npre=NPRE, n_gh=24):
    """TEB de la QPSK OFDM, porteuse par porteuse, quand le canal est estimé par la moyenne de npre symboles
    d'apprentissage. Après la FFT (divisée par N), Y = H X + n avec n ~ CN(0, 1/(2 Eb/N0)) ; l'estimation vaut
    H + e avec e ~ CN(0, 1/(2 npre Eb/N0)). Pour une erreur e donnée, chaque bit de la QPSK se trompe avec la
    probabilité Q(a Re(H X conj(H + e)) / sqrt(|H + e|^2 / (4 Eb/N0))) ; la moyenne sur e est calculée
    exactement par une quadrature de Gauss-Hermite à deux dimensions (n_gh x n_gh points)."""
    t, w = np.polynomial.hermite.hermgauss(n_gh)
    u, v = np.meshgrid(np.sqrt(2) * t, np.sqrt(2) * t); W = np.outer(w, w) / np.pi
    se = np.sqrt(1 / (2 * npre * ebn0_lin) / 2)
    e = se * (u + 1j * v); X = (1 + 1j) / np.sqrt(2); p = 0.0
    for h in H(taps):
        Hh = h + e; s = h * X * np.conj(Hh); sd = np.abs(Hh) / (2 * np.sqrt(ebn0_lin))
        p += np.sum(W * 0.5 * (Q(s.real / sd) + Q(s.imag / sd)))
    return float(p / len(OCC))


def teb_canal_connu(taps, ebn0_lin):
    """TEB de la QPSK OFDM avec un canal parfaitement connu"""
    return float(np.mean(Q(np.sqrt(2 * np.abs(H(taps)) ** 2 * ebn0_lin))))


def table(pts):
    return '[' + ', '.join(f'complex({p.real!r}, {p.imag!r})' for p in pts) + ']'
