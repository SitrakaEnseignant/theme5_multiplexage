"""
Bloc Python intégré « egaliseur » : estimation du canal sur les symboles d'apprentissage,
égalisation à un coefficient par porteuse et extraction des porteuses de données.
Entrée : symboles OFDM après la FFT (vecteurs de fft_len valeurs, porteuse 0 au centre).
Sortie 0 : les nd symboles de données égalisés de chaque symbole OFDM de données.
Sortie 1 : l'estimation du canal sur les fft_len porteuses, répétée pour chaque symbole de données.
Chaque trame comprend npre symboles d'apprentissage identiques puis nsym symboles de données.
"""
import numpy as np
from gnuradio import gr


class blk(gr.basic_block):
    def __init__(self, fft_len=64, npre=2, nsym=20, preambule=[0.0] * 64, porteuses=[0] * 48):
        self.nd = len(porteuses)
        gr.basic_block.__init__(self, name='egaliseur', in_sig=[(np.complex64, fft_len)],
                                out_sig=[(np.complex64, self.nd), (np.complex64, fft_len)])
        self.npre, self.nsym = npre, nsym
        self.pre = np.array(preambule, dtype=complex)
        self.occ = np.array(porteuses, dtype=int) + fft_len // 2
        self.k, self.acc, self.h = 0, np.zeros(fft_len, complex), np.ones(fft_len, complex)

    def general_work(self, input_items, output_items):
        x, y0, y1 = input_items[0], output_items[0], output_items[1]
        i = o = 0
        while i < len(x) and o < len(y0):
            if self.k < self.npre:                       # symbole d'apprentissage : accumulation
                self.acc = self.acc + x[i]
                if self.k == self.npre - 1:              # estimation H = moyenne(Y) / P sur les porteuses utilisées
                    p = self.pre != 0
                    self.h = np.ones(len(self.acc), complex)
                    self.h[p] = self.acc[p] / self.npre / self.pre[p]
            else:                                        # symbole de données : égalisation Z = Y / H
                y0[o] = (x[i] / self.h)[self.occ]
                y1[o] = self.h
                o += 1
            self.k += 1
            if self.k == self.npre + self.nsym:
                self.k, self.acc = 0, np.zeros(len(self.acc), complex)
            i += 1
        self.consume(0, i)
        return o
