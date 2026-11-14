# theme5_multiplexage

Bloc 2, Thème 5 : multiplexage. Le multiplexage temporel (TDM) et le multiplexage fréquentiel (FDM) sont
étudiés dans le cours par leurs schémas synoptiques ; le multiplexage par répartition orthogonale de la
fréquence (OFDM) est simulé dans MATLAB Simulink et dans GNU Radio 3.10, en simulation pure (aucun matériel
de radio logicielle).

| Paramètre | Valeur |
|---|---|
| Cadence d'échantillonnage | 1 MHz, comme la leçon 30 du cours Learn SDR |
| FFT | 64 porteuses, écart de 15 625 Hz |
| Préfixe cyclique | 16 échantillons (16 µs) |
| Porteuses | disposition du 802.11a : 48 porteuses de données en QPSK, 4 pilotes, porteuse 0 vide |
| Trame | 2 symboles d'apprentissage identiques puis 20 symboles de données |
| Récepteur | retrait du préfixe, FFT, canal estimé par la moyenne des 2 symboles d'apprentissage, égalisation Y/H |
| Canaux | bruit seul ; écho de 0,9 à 3 µs (plus court que le préfixe) ; écho de 0,9 à 20 µs (plus long que le préfixe) |

## Contenu

- `En_ligne/ofdm.py` : paramètres communs et formules théoriques (TEB avec canal estimé).
- `outils/generer_grc.py` : génère `GNU_Radio/ofdm_tx_rx_sim.grc` et `Simulink/parametres_ofdm.json`.
- `outils/egaliseur_epy.py` : code du bloc Python intégré « egaliseur » du graphe de flux.
- `Simulink/construire_et_verifier_modeles.m` : construit `mod_ofdm`, l'exécute et vérifie les résultats.
- `En_ligne/executer_grc_ci.py` : exécution sans écran du graphe de flux, mesures et figures.
- `.github/workflows/simulations.yml` : exécution automatique ; résultats recopiés dans `Simulink/modeles`,
  `Simulink/resultats`, `GNU_Radio/resultats` et dans la branche `resultats`.
