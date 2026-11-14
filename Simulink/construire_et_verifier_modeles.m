%% CONSTRUIRE_ET_VERIFIER_MODELES - Bloc 2, Thème 5 : multiplexage OFDM
% Construit le modèle Simulink mod_ofdm (simulation pure, sans matériel), qui reprend étage par étage le
% graphe de flux GNU Radio GNU_Radio/ofdm_tx_rx_sim.grc : source et table QPSK, allocateur de porteuses
% (préambule, données, pilotes), IFFT, préfixe cyclique, canal à trajets multiples, bruit gaussien complexe,
% retrait du préfixe, FFT, estimation du canal et égalisation, décision. Le modèle fonctionne à la cadence
% des symboles OFDM : chaque pas de calcul traite un symbole complet (vecteurs de 48, 64 ou 80 valeurs).
% Chaque élément des vecteurs aléatoires a sa propre graine, sinon les 48 (ou 80) valeurs d'un pas seraient identiques.
% Le script exécute le modèle sans bruit pour trois canaux, puis pour plusieurs Eb/N0, et vérifie le
% spectre, la puissance, le canal estimé et le taux d'erreur binaire comparé à la théorie.
% Les paramètres sont lus dans parametres_ofdm.json, écrit par outils/generer_grc.py.
clear; clc;
ici = fileparts(mfilename('fullpath')); cd(ici);
if exist(fullfile(ici, 'journal_matlab.txt'), 'file'), delete(fullfile(ici, 'journal_matlab.txt')); end
diary(fullfile(ici, 'journal_matlab.txt'));
fprintf('MATLAB %s\n', version);
P = jsondecode(fileread(fullfile(ici, 'parametres_ofdm.json')));
P.OCC = P.OCC(:).'; P.PIL = P.PIL(:).'; P.PILS = P.PILS(:).'; P.USED = P.USED(:).'; P.PRE = P.PRE(:).';
P.NS = P.N + P.CP; P.ND = numel(P.OCC); P.Tsym = P.NS / P.fs; P.pts = complex(P.QPSK_re(:).', P.QPSK_im(:).'); P.LAB = P.LAB(:).';
P.DF = P.fs / P.N;
canaux = {'bbag', 'echo_court', 'echo_long'};
EBN0 = struct('bbag', [0 2 4 6 8 10], 'echo_court', [0 4 8 12 16 20]);
setappdata(0, 'liste_blocs', {});
mdl = construire(P);
fid = fopen(fullfile(ici, [mdl '_blocs.json']), 'w'); fprintf(fid, '%s', jsonencode(getappdata(0, 'liste_blocs'), 'PrettyPrint', true)); fclose(fid);
save_system(mdl, fullfile(ici, [mdl '.slx']));
try, print(['-s' mdl], '-dpng', '-r150', fullfile(ici, [mdl '_schema.png'])); catch err, fprintf('Schema non exporte : %s\n', err.message); end
R = struct(); echecs = {}; D = struct();
% 1. sans bruit, trois canaux
for i = 1:numel(canaux)
    c = canaux{i}; h = P.canaux.(c)(:).';
    so = simuler(mdl, P, 10, h, 0); [es, ns, eb, nb, e, z, ref] = mesurer(so, P);
    he = lirev(so, 'hest', P.N); he = he(end, :) / P.N; ht = reponse(h, P);
    eh = max(abs(he(P.OCC + 33) - ht(P.OCC + 33)));
    R.(c).sans_bruit = struct('erreurs_symboles', es, 'symboles', ns, 'TEB', eb / nb, 'evm', e, 'ecart_max_canal_estime', eh);
    fprintf('%s sans bruit : %d erreurs sur %d symboles, EVM %.2e, ecart canal %.2e\n', c, es, ns, e, eh);
    if strcmp(c, 'echo_long')
        if e < 0.1, echecs{end + 1} = sprintf('echo_long : EVM %.3f alors qu''une degradation est attendue', e); end %#ok<AGROW>
    elseif es > 0 || e > 1e-3 || eh > 1e-3
        echecs{end + 1} = sprintf('%s : %d erreurs, EVM %.2e, ecart canal %.2e sans bruit', c, es, e, eh); %#ok<AGROW>
    end
    D.(c) = so;
end
% spectre, puissance moyenne et PAPR du signal émis (canal idéal, sans bruit)
s = lirev(D.bbag, 's', P.NS); s = reshape(s.', 1, []);
[f, Pw] = densite(s, P.fs); frac = sum(Pw(abs(f) <= 26.5 * P.DF)) / sum(Pw); pm = mean(abs(s) .^ 2);
X = reshape(s(1:floor(numel(s) / P.NS) * P.NS), P.NS, []).'; X = X(:, P.CP + 1:end);
papr = 10 * log10(max(abs(X) .^ 2, [], 2) / pm); p8 = mean(papr > 8); p8t = 1 - (1 - exp(-10 ^ 0.8)) ^ 52;
R.signal = struct('part_puissance_bande', frac, 'puissance_moyenne', pm, 'puissance_attendue', 52, ...
    'papr_moyen_dB', mean(papr), 'P_papr_sup_8dB', p8, 'P_papr_sup_8dB_approx', p8t);
fprintf('Spectre : %.4f de la puissance dans la bande ; puissance moyenne %.2f (attendue 52) ; P(PAPR > 8 dB) = %.3f (approx. %.3f)\n', frac, pm, p8, p8t);
if frac < 0.98, echecs{end + 1} = sprintf('spectre : %.3f de la puissance dans la bande', frac); end
if abs(pm / 52 - 1) > 0.05, echecs{end + 1} = sprintf('puissance moyenne %.2f au lieu de 52', pm); end
% 2. taux d'erreur binaire
for c = {'bbag', 'echo_court'}
    c = c{1}; h = P.canaux.(c)(:).'; T = struct('EbN0_dB', {}, 'TEB', {}, 'erreurs', {}, 'bits', {}, 'theorie', {}, 'canal_connu', {}, 'conforme', {});
    for e_db = EBN0.(c)
        g = 10 ^ (e_db / 10); s2 = P.N / (2 * g);
        so = simuler(mdl, P, 100, h, s2); [~, ~, eb, nb] = mesurer(so, P);
        th = teb_theorique(h, g, P); mes = eb / nb; ok = eb < 50 || (mes / th > 1 / 1.5 && mes / th < 1.5);
        T(end + 1) = struct('EbN0_dB', e_db, 'TEB', mes, 'erreurs', eb, 'bits', nb, 'theorie', th, ...
            'canal_connu', teb_canal_connu(h, g, P), 'conforme', ok); %#ok<AGROW>
        fprintf('%s Eb/N0 = %2d dB : TEB %.3e (%d erreurs), theorie %.3e  %s\n', c, e_db, mes, eb, th, string(ok));
        if ~ok, echecs{end + 1} = sprintf('%s : Eb/N0 %d dB, TEB %.3g, theorie %.3g', c, e_db, mes, th); end %#ok<AGROW>
        if strcmp(c, 'echo_court') && e_db == 20, D.echo20 = so; end
    end
    R.(c).taux = T;
end
figures(P, D, R, s, f, Pw, papr, ici);
fid = fopen(fullfile(ici, 'mesures_simulink.json'), 'w'); fprintf(fid, '%s', jsonencode(R, 'PrettyPrint', true)); fclose(fid);
close_system(mdl, 0);
fprintf('\nControles en echec : %d\n', numel(echecs));
for k = 1:numel(echecs), fprintf('  ECHEC %s\n', echecs{k}); end
diary off;
if ~isempty(echecs), error('Des controles ont echoue (voir journal_matlab.txt).'); end

%% ------------------------------------------------------------------ construction du modèle
function b = ajouter(mdl, lib, nom, col, lig, varargin)
x = 30 + 150 * col; y = 30 + 100 * lig;
b = [mdl '/' nom];
L = getappdata(0, 'liste_blocs'); L{end + 1} = struct('bloc', nom, 'bibliotheque', strrep(lib, 'simulink/', ''), 'parametres', {varargin}); setappdata(0, 'liste_blocs', L);
add_block(lib, b, 'Position', [x y x + 90 y + 50], varargin{:});
end

function relier(mdl, a, b)
add_line(mdl, a, b, 'autorouting', 'on');
end

function fonction(mdl, nom, col, lig, code)
ajouter(mdl, 'simulink/User-Defined Functions/MATLAB Function', nom, col, lig);
ch = sfroot().find('-isa', 'Stateflow.EMChart', 'Path', [mdl '/' nom]); ch.Script = code;
end

function v = vec(x)
% écriture d'un vecteur ligne pour le code des blocs MATLAB Function
v = mat2str(x(:).', 17);
end

function mdl = construire(P)
mdl = 'mod_ofdm'; B = 'simulink/';
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl);
set_param(mdl, 'SolverType', 'Fixed-step', 'Solver', 'FixedStepDiscrete', 'FixedStep', num2str(P.Tsym, 17), 'StopTime', '1');
hws = get_param(mdl, 'ModelWorkspace');
v = {'Tsym', P.Tsym; 'sigma2', 0; 'taps', 1; 'graine', 7; 'ptsI', real(P.pts); 'ptsQ', imag(P.pts); 'ND', P.ND; 'NS', P.NS};
for k = 1:size(v, 1), hws.assignin(v{k, 1}, v{k, 2}); end
% émetteur
ajouter(mdl, [B 'Sources/Uniform Random Number'], 'src', 0, 1, 'Minimum', 'zeros(1, ND)', 'Maximum', '4 * ones(1, ND)', 'Seed', 'graine + (1:ND)', 'SampleTime', 'Tsym');
ajouter(mdl, [B 'Math Operations/Rounding Function'], 'indice', 1, 1, 'Operator', 'floor');
ajouter(mdl, [B 'Lookup Tables/Direct Lookup Table (n-D)'], 'table_I', 2, 0.5, 'NumberOfTableDimensions', '1', 'Table', 'ptsI');
ajouter(mdl, [B 'Lookup Tables/Direct Lookup Table (n-D)'], 'table_Q', 2, 1.5, 'NumberOfTableDimensions', '1', 'Table', 'ptsQ');
ajouter(mdl, [B 'Math Operations/Real-Imag to Complex'], 'symboles_tx', 3, 1, 'Input', 'Real and imag');
ajouter(mdl, [B 'Sinks/To Workspace'], 'sym', 2, 2.6, 'VariableName', 'sym', 'SaveFormat', 'Array');
relier(mdl, 'src/1', 'indice/1'); relier(mdl, 'indice/1', 'table_I/1'); relier(mdl, 'indice/1', 'table_Q/1'); relier(mdl, 'indice/1', 'sym/1');
relier(mdl, 'table_I/1', 'symboles_tx/1'); relier(mdl, 'table_Q/1', 'symboles_tx/2');
fonction(mdl, 'allocateur', 4, 1, code_allocateur(P));
fonction(mdl, 'ifft_tx', 5, 1, sprintf(['function x = ifft_ofdm(X)\n%% IFFT non normalisee, porteuse 0 au centre en entree (comme le bloc FFT inverse de GNU Radio)\n' ...
    'x = %d * ifft(ifftshift(X(:).''));\nend\n'], P.N));
fonction(mdl, 'prefixe', 6, 1, sprintf(['function y = prefixe(x)\n%% prefixe cyclique : copie des %d derniers echantillons en tete du symbole\n' ...
    'x = x(:).''; y = [x(%d:%d), x];\nend\n'], P.CP, P.N - P.CP + 1, P.N));
ajouter(mdl, [B 'Sinks/To Workspace'], 's', 7, 0, 'VariableName', 's', 'SaveFormat', 'Array');
ajouter(mdl, [B 'Sinks/To Workspace'], 'donnees', 5, 2.4, 'VariableName', 'donnees', 'SaveFormat', 'Array');
relier(mdl, 'symboles_tx/1', 'allocateur/1'); relier(mdl, 'allocateur/1', 'ifft_tx/1'); relier(mdl, 'allocateur/2', 'donnees/1');
relier(mdl, 'ifft_tx/1', 'prefixe/1'); relier(mdl, 'prefixe/1', 's/1');
% canal
ajouter(mdl, [B 'Sources/Constant'], 'taps', 6, 2.2, 'Value', 'taps');
fonction(mdl, 'canal_trajets', 7, 1, code_canal_trajets());
ajouter(mdl, [B 'Sources/Random Number'], 'bruit_I', 7, 2.4, 'Mean', 'zeros(1, NS)', 'Variance', 'sigma2 / 2', 'Seed', '1000 + (1:NS)', 'SampleTime', 'Tsym');
ajouter(mdl, [B 'Sources/Random Number'], 'bruit_Q', 7, 3.4, 'Mean', 'zeros(1, NS)', 'Variance', 'sigma2 / 2', 'Seed', '2000 + (1:NS)', 'SampleTime', 'Tsym');
fonction(mdl, 'canal', 8, 1.6, sprintf('function r = canal(x, nI, nQ)\n%% addition du bruit gaussien complexe\nr = x(:).'' + complex(nI(:).'', nQ(:).'');\nend\n'));
relier(mdl, 'prefixe/1', 'canal_trajets/1'); relier(mdl, 'taps/1', 'canal_trajets/2'); relier(mdl, 'canal_trajets/1', 'canal/1');
relier(mdl, 'bruit_I/1', 'canal/2'); relier(mdl, 'bruit_Q/1', 'canal/3');
% récepteur
fonction(mdl, 'retrait_prefixe', 9, 1.6, sprintf('function y = retrait_prefixe(r)\n%% retrait du prefixe cyclique (equivalent du bloc Keep M in N)\nr = r(:).''; y = r(%d:%d);\nend\n', P.CP + 1, P.NS));
fonction(mdl, 'fft_rx', 10, 1.6, sprintf('function Y = fft_ofdm(y)\n%% FFT non normalisee, porteuse 0 au centre en sortie\nY = fftshift(fft(y(:).''));\nend\n'));
fonction(mdl, 'egaliseur', 11, 1.6, code_egaliseur(P));
fonction(mdl, 'decision', 12, 1.6, code_decision(P));
ajouter(mdl, [B 'Sinks/To Workspace'], 'Ybrut', 11, 3.2, 'VariableName', 'Ybrut', 'SaveFormat', 'Array');
ajouter(mdl, [B 'Sinks/To Workspace'], 'z', 13, 0.6, 'VariableName', 'z', 'SaveFormat', 'Array');
ajouter(mdl, [B 'Sinks/To Workspace'], 'hest', 13, 1.6, 'VariableName', 'hest', 'SaveFormat', 'Array');
ajouter(mdl, [B 'Sinks/To Workspace'], 'valide', 13, 2.6, 'VariableName', 'valide', 'SaveFormat', 'Array');
ajouter(mdl, [B 'Sinks/To Workspace'], 'dec', 13, 3.6, 'VariableName', 'dec', 'SaveFormat', 'Array');
relier(mdl, 'canal/1', 'retrait_prefixe/1'); relier(mdl, 'retrait_prefixe/1', 'fft_rx/1'); relier(mdl, 'fft_rx/1', 'egaliseur/1');
relier(mdl, 'fft_rx/1', 'Ybrut/1'); relier(mdl, 'egaliseur/1', 'decision/1'); relier(mdl, 'egaliseur/1', 'z/1');
relier(mdl, 'egaliseur/2', 'hest/1'); relier(mdl, 'egaliseur/3', 'valide/1'); relier(mdl, 'decision/1', 'dec/1');
end

function c = code_allocateur(P)
c = sprintf(['function [X, donnees] = allocateur(d)\n' ...
    '%% repartition sur les porteuses (porteuse 0 a l''indice 33) : %d symboles d''apprentissage puis %d symboles\n' ...
    '%% de donnees par trame, %d porteuses de donnees et %d pilotes (equivalent du bloc OFDM Carrier Allocator)\n' ...
    'persistent k\nif isempty(k), k = 0; end\n' ...
    'pre = %s;\nocc = %s;\npil = %s;\npils = %s;\n' ...
    'X = complex(zeros(1, %d));\n' ...
    'if k < %d\n    X = complex(pre); donnees = 0;\nelse\n    d = d(:).''; X(occ) = d; X(pil) = pils; donnees = 1;\nend\n' ...
    'k = k + 1; if k == %d, k = 0; end\nend\n'], P.NPRE, P.NSYM, P.ND, numel(P.PIL), vec(P.PRE), vec(P.OCC + 33), vec(P.PIL + 33), vec(P.PILS), ...
    P.N, P.NPRE, P.NPRE + P.NSYM);
end

function c = code_canal_trajets()
c = sprintf(['function y = canal_trajets(x, h)\n' ...
    '%% canal a trajets multiples : filtre a reponse impulsionnelle finie de coefficients h, avec memoire\n' ...
    '%% d''un symbole OFDM au suivant (equivalent du bloc FIR Filter)\n' ...
    'persistent mem\nh = h(:).''; x = x(:).''; L = numel(h);\n' ...
    'if isempty(mem), mem = complex(zeros(1, max(L - 1, 1))); end\n' ...
    'u = [mem(1:L - 1), x]; y = complex(zeros(1, numel(x)));\n' ...
    'for n = 1:numel(x)\n    acc = complex(0);\n    for k = 1:L\n        acc = acc + h(k) * u(n + L - k);\n    end\n    y(n) = acc;\nend\n' ...
    'if L > 1, mem(1:L - 1) = u(end - L + 2:end); end\nend\n']);
end

function c = code_egaliseur(P)
c = sprintf(['function [z, hest, valide] = egaliseur(Y)\n' ...
    '%% estimation du canal par la moyenne des %d symboles d''apprentissage, puis egalisation Z = Y / H\n' ...
    '%% sur chaque porteuse et extraction des %d porteuses de donnees (meme code que le bloc Python de GNU Radio)\n' ...
    'persistent k acc h\n' ...
    'if isempty(k), k = 0; acc = complex(zeros(1, %d)); h = complex(ones(1, %d)); end\n' ...
    'pre = %s;\nused = %s;\nocc = %s;\n' ...
    'Y = Y(:).''; z = complex(zeros(1, %d)); valide = 0;\n' ...
    'if k < %d\n    acc = acc + Y;\n    if k == %d\n        h = complex(ones(1, %d)); h(used) = acc(used) / %d ./ pre(used);\n    end\n' ...
    'else\n    zz = Y ./ h; z = zz(occ); valide = 1;\nend\n' ...
    'hest = h;\nk = k + 1; if k == %d, k = 0; acc = complex(zeros(1, %d)); end\nend\n'], ...
    P.NPRE, P.ND, P.N, P.N, vec(P.PRE), vec(P.USED + 33), vec(P.OCC + 33), P.ND, P.NPRE, P.NPRE - 1, P.N, P.NPRE, P.NPRE + P.NSYM, P.N);
end

function c = code_decision(P)
c = sprintf(['function idx = decision(z)\n%% decision au plus proche point de la constellation QPSK (Constellation Decoder)\n' ...
    'pts = complex(%s, %s);\nz = z(:).''; idx = zeros(1, numel(z));\n' ...
    'for i = 1:numel(z)\n    [~, m] = min(abs(z(i) - pts)); idx(i) = m - 1;\nend\nend\n'], vec(real(P.pts)), vec(imag(P.pts)));
end

%% ------------------------------------------------------------------ exécution et mesures
function so = simuler(mdl, P, ntrames, h, s2)
hws = get_param(mdl, 'ModelWorkspace'); hws.assignin('taps', h); hws.assignin('sigma2', s2);
n = ntrames * (P.NPRE + P.NSYM);
so = sim(mdl, 'ReturnWorkspaceOutputs', 'on', 'StopTime', num2str((n - 1) * P.Tsym, 17));
end

function v = lirev(so, nom, n)
% signal enregistré par To Workspace, ramené à une matrice (pas de calcul x n)
v = so.get(nom);
if ndims(v) == 3, v = reshape(v, size(v, 1) * size(v, 2), size(v, 3)).'; end
if size(v, 2) ~= n && size(v, 1) == n, v = v.'; end
end

function [es, ns, eb, nb, e, z, ref] = mesurer(so, P)
sym = lirev(so, 'sym', P.ND); don = lirev(so, 'donnees', 1); dec = lirev(so, 'dec', P.ND);
val = lirev(so, 'valide', 1); z = lirev(so, 'z', P.ND);
a = round(sym(don == 1, :)); b = round(dec(val == 1, :)); z = z(val == 1, :);
n = min(size(a, 1), size(b, 1)); a = a(1:n, :); b = b(1:n, :); z = z(1:n, :);
es = sum(a(:) ~= b(:)); ns = numel(a); x = bitxor(P.LAB(a(:) + 1), P.LAB(b(:) + 1));
eb = sum(x == 1 | x == 2) + 2 * sum(x == 3); nb = 2 * ns;
ref = P.pts(a + 1); e = sqrt(mean(abs(z(:) - ref(:)) .^ 2));
end

function H = reponse(h, P)
% réponse en fréquence du canal sur les 64 porteuses, porteuse 0 à l'indice 33
k = -P.N / 2:P.N / 2 - 1; H = zeros(1, P.N);
for i = 1:P.N, H(i) = sum(h .* exp(-2j * pi * k(i) * (0:numel(h) - 1) / P.N)); end
end

function q = Qf(x), q = 0.5 * erfc(x / sqrt(2)); end

function p = teb_canal_connu(h, g, P)
H = reponse(h, P); H = H(P.OCC + 33); p = mean(Qf(sqrt(2 * abs(H) .^ 2 * g)));
end

function p = teb_theorique(h, g, P)
% TEB de la QPSK OFDM avec un canal estimé par la moyenne de NPRE symboles d'apprentissage :
% moyenne exacte sur l'erreur d'estimation, par quadrature de Gauss-Hermite à deux dimensions (24 x 24)
n = 24; i = 1:n - 1; J = diag(sqrt(i / 2), 1) + diag(sqrt(i / 2), -1); [V, Dg] = eig(J);
t = diag(Dg); w = sqrt(pi) * V(1, :).' .^ 2;
[u, v] = meshgrid(sqrt(2) * t, sqrt(2) * t); W = (w * w.') / pi;
se = sqrt(1 / (2 * P.NPRE * g) / 2); e = se * (u + 1j * v); X = (1 + 1j) / sqrt(2);
H = reponse(h, P); H = H(P.OCC + 33); p = 0;
for k = 1:numel(H)
    Hh = H(k) + e; s = H(k) * X * conj(Hh); sd = abs(Hh) / (2 * sqrt(g));
    p = p + sum(sum(W .* 0.5 .* (Qf(real(s) ./ sd) + Qf(imag(s) ./ sd))));
end
p = p / numel(H);
end

function [f, Pw] = densite(s, fs)
Nf = 1024; w = 0.5 - 0.5 * cos(2 * pi * (0:Nf - 1) / Nf); Pw = zeros(1, Nf); k = 0;
for i = 1:Nf / 2:numel(s) - Nf
    X = fftshift(fft(s(i:i + Nf - 1) .* w)); Pw = Pw + abs(X) .^ 2; k = k + 1;
end
Pw = Pw / k; f = (-Nf / 2:Nf / 2 - 1) * fs / Nf;
end

function figures(P, D, R, s, f, Pw, papr, ici)
fh = figure('Visible', 'off', 'Position', [0 0 1700 1000]);
subplot(2, 3, 1); t = (0:3 * P.NS - 1) / P.fs * 1e6; s0 = s(2 * P.NS + 1:5 * P.NS);
hold on; for k = 0:2, patch([k * P.NS k * P.NS + P.CP k * P.NS + P.CP k * P.NS], [-30 -30 30 30], [1 0.8 0.5], 'EdgeColor', 'none'); end
h1 = plot(t, real(s0)); h2 = plot(t, abs(s0)); grid on; xlabel('temps (\mus)'); ylim([-30 30]);
title('signal emis ; prefixes cycliques en orange'); legend([h1 h2], {'partie reelle', 'module'}, 'Location', 'southwest');
subplot(2, 3, 2); plot(f / 1e3, 10 * log10(Pw / max(Pw) + 1e-12)); ylim([-60 3]); grid on; xlabel('frequence (kHz)'); ylabel('dB'); title('spectre du signal emis');
subplot(2, 3, 3); zz = sort(papr); semilogy(zz, 1 - (0:numel(zz) - 1) / numel(zz), '.'); hold on;
zt = linspace(4, 12, 50); semilogy(zt, 1 - (1 - exp(-10 .^ (zt / 10))) .^ 52); ylim([1e-3 1.2]); grid on;
xlabel('PAPR (dB)'); ylabel('P(PAPR > z)'); title('rapport puissance crete sur puissance moyenne'); legend('mesure', 'approximation', 'Location', 'southwest');
so = D.echo20; Yb = lirev(so, 'Ybrut', P.N); val = lirev(so, 'valide', 1); Yb = Yb(val == 1, P.OCC + 33) / P.N;
subplot(2, 3, 4); plot(real(Yb(:)), imag(Yb(:)), '.', 'MarkerSize', 3); axis equal; grid on; title('porteuses avant egalisation, echo court, 20 dB');
z = lirev(so, 'z', P.ND); z = z(val == 1, :);
subplot(2, 3, 5); plot(real(z(:)), imag(z(:)), '.', 'MarkerSize', 3); axis equal; xlim([-2 2]); ylim([-2 2]); grid on; title('apres egalisation, echo court, 20 dB');
subplot(2, 3, 6); hold on; k = -P.N / 2:P.N / 2 - 1; co = {[0 0.45 0.74], [0.85 0.33 0.1]}; nm = {'echo_court', 'echo_long'};
for i = 1:2
    he = lirev(D.(nm{i}), 'hest', P.N); he = he(end, :) / P.N; ht = reponse(P.canaux.(nm{i})(:).', P);
    plot(k(P.OCC + 33), abs(he(P.OCC + 33)), 'o', 'Color', co{i}, 'MarkerSize', 4); plot(k, abs(ht), '-', 'Color', co{i});
end
grid on; xlabel('indice de porteuse'); title('module du canal estime, sans bruit');
legend('echo court estime', 'echo court theorique', 'echo long estime', 'echo long theorique', 'Location', 'south');
sgtitle(['mod\_ofdm  -  MATLAB R' version('-release')]);
exportgraphics(fh, fullfile(ici, 'mod_ofdm_resultat.png'), 'Resolution', 110); close(fh);
fh = figure('Visible', 'off', 'Position', [0 0 1000 700]); co = lines(2); nm = {'bbag', 'echo_court'}; hold on;
for i = 1:2
    T = R.(nm{i}).taux; e = [T.EbN0_dB]; y = [T.TEB]; y(y == 0) = NaN; ee = linspace(min(e), max(e), 40); h = P.canaux.(nm{i})(:).';
    semilogy(ee, arrayfun(@(x) teb_theorique(h, 10 ^ (x / 10), P), ee), '-', 'Color', co(i, :));
    semilogy(ee, arrayfun(@(x) teb_canal_connu(h, 10 ^ (x / 10), P), ee), '--', 'Color', co(i, :));
    semilogy(e, y, 'o', 'Color', co(i, :));
end
set(gca, 'YScale', 'log'); ylim([1e-5 0.5]); grid on; xlabel('Eb/N0 (dB)'); ylabel('TEB');
legend('bbag : theorie, canal estime', 'bbag : canal connu', 'bbag : mesure Simulink', 'echo court : theorie, canal estime', ...
    'echo court : canal connu', 'echo court : mesure Simulink', 'Location', 'southwest');
title(['TEB de la QPSK OFDM mesure dans Simulink et theorie - MATLAB R' version('-release')]);
exportgraphics(fh, fullfile(ici, 'teb_simulink.png'), 'Resolution', 110); close(fh);
end
