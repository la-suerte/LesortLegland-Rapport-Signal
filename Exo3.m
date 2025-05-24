%% Exo - Filtrage Audio avec Filtres IIR
% Objectif: Éliminer les bruits parasites d'un enregistrement musical
% Méthode: Application de filtres coupe-bande en cascade
% Auteur: Étudiant - Traitement du Signal

% Initialisation de l'environnement de travail
clear variables; clc; close all;

%% Partie 1: Chargement et analyse du fichier audio
fprintf('=== ÉTAPE 1: CHARGEMENT DES DONNÉES ===\n');

% Lecture du fichier audio Mozart avec bruit
nom_fichier = 'Mozart_Bruit.wav';
[signal_mozart, freq_echantillonnage] = audioread(nom_fichier);

% Conversion stéréo vers mono si nécessaire
if size(signal_mozart, 2) > 1
    signal_mozart = signal_mozart(:,1);  % Prendre seulement le canal gauche
    fprintf('Conversion stéréo -> mono effectuée\n');
end

% Récupération des métadonnées du fichier
infos_fichier = audioinfo(nom_fichier);
nb_echantillons = length(signal_mozart);
duree_enregistrement = nb_echantillons / freq_echantillonnage;
pas_temporel = 1 / freq_echantillonnage;
bits_quantification = infos_fichier.BitsPerSample;

% Création de l'axe temporel
temps_signal = (0:nb_echantillons-1) * pas_temporel;

% Affichage des caractéristiques du signal
fprintf('Caractéristiques du signal audio:\n');
fprintf('- Période d''échantillonnage: %.6f s\n', pas_temporel);
fprintf('- Fréquence d''échantillonnage: %d Hz\n', freq_echantillonnage);
fprintf('- Résolution: %d bits\n', bits_quantification);
fprintf('- Durée totale: %.2f secondes\n', duree_enregistrement);
fprintf('- Nombre d''échantillons: %d\n', nb_echantillons);

%% Partie 2: Analyse spectrale du signal original
fprintf('\n=== ÉTAPE 2: ANALYSE SPECTRALE ===\n');

% Calcul de la FFT du signal d'entrée
transformation_fourier = fft(signal_mozart);
magnitude_spectre = abs(transformation_fourier) / nb_echantillons;

% Construction du spectre unilatéral (fréquences positives)
spectre_positif = magnitude_spectre(1:floor(nb_echantillons/2)+1);
spectre_positif(2:end-1) = 2 * spectre_positif(2:end-1);

% Vecteur des fréquences correspondant
resolution_spectrale = freq_echantillonnage / nb_echantillons;
vecteur_freq = (0:floor(nb_echantillons/2)) * resolution_spectrale;

fprintf('Résolution spectrale: %.3f Hz\n', resolution_spectrale);

%% Partie 3: Conception et application des filtres
fprintf('\n=== ÉTAPE 3: FILTRAGE DES PARASITES ===\n');

% Fréquence de Nyquist (limite théorique)
freq_nyquist = freq_echantillonnage / 2;

% Liste des fréquences parasites identifiées lors de l'analyse préalable
parasites_freq = [0.03, 117.01, 160.93, 333.19, 1397.01, 5658.01];
largeur_bande = 25;  % Largeur de la bande coupée (±25 Hz)

% Initialisation du signal de sortie
signal_nettoye = signal_mozart;

% Coefficients du filtre global (pour analyse finale)
numerateur_global = 1;
denominateur_global = 1;

% Application successive des filtres coupe-bande
for num_filtre = 1:length(parasites_freq)
    freq_centrale = parasites_freq(num_filtre);
    
    % Ignorer les fréquences trop basses (problèmes numériques)
    if freq_centrale < 1
        fprintf('Fréquence %.2f Hz ignorée (trop basse)\n', freq_centrale);
        continue;
    end
    
    % Calcul des bornes de la bande à éliminer
    freq_inf = max(0.001, (freq_centrale - largeur_bande) / freq_nyquist);
    freq_sup = min(0.999, (freq_centrale + largeur_bande) / freq_nyquist);
    
    % Conception du filtre Butterworth coupe-bande d'ordre 2
    [coef_num, coef_den] = butter(2, [freq_inf, freq_sup], 'stop');
    
    % Application du filtre au signal
    signal_nettoye = filter(coef_num, coef_den, signal_nettoye);
    
    % Mise à jour des coefficients globaux (produit de convolution)
    numerateur_global = conv(numerateur_global, coef_num);
    denominateur_global = conv(denominateur_global, coef_den);
    
    % Information sur le filtre appliqué
    fprintf('Filtre n°%d appliqué:\n', num_filtre);
    fprintf('  - Fréquence centrale: %.2f Hz\n', freq_centrale);
    fprintf('  - Bande supprimée: [%.2f - %.2f] Hz\n', ...
            freq_centrale-largeur_bande, freq_centrale+largeur_bande);
end

%% Partie 4: Post-traitement du signal filtré
fprintf('\n=== ÉTAPE 4: NORMALISATION ===\n');

% Normalisation pour éviter la saturation
amplitude_max = max(abs(signal_nettoye));
signal_final = signal_nettoye / amplitude_max;

fprintf('Facteur de normalisation appliqué: %.4f\n', 1/amplitude_max);

% Calcul du spectre du signal filtré
spectre_filtre_fft = abs(fft(signal_final)) / nb_echantillons;
spectre_filtre_final = spectre_filtre_fft(1:floor(nb_echantillons/2)+1);
spectre_filtre_final(2:end-1) = 2 * spectre_filtre_final(2:end-1);

%% Partie 5: Visualisation des résultats
fprintf('\n=== ÉTAPE 5: GÉNÉRATION DES GRAPHIQUES ===\n');

% Graphique 1: Comparaison temporelle
figure('Name', 'Analyse temporelle', 'Position', [100 100 800 600]);

subplot(2, 1, 1);
plot(temps_signal, signal_mozart, 'Color', [0.2 0.4 0.8], 'LineWidth', 0.7);
title('Signal Audio Original (avec parasites)');
xlabel('Temps (secondes)');
ylabel('Amplitude');
grid on; grid minor;

subplot(2, 1, 2);
plot(temps_signal, signal_final, 'Color', [0.8 0.2 0.2], 'LineWidth', 1);
title('Signal Audio Nettoyé (après filtrage)');
xlabel('Temps (secondes)');
ylabel('Amplitude normalisée');
grid on; grid minor;

% Graphique 2: Comparaison spectrale
figure('Name', 'Analyse spectrale', 'Position', [200 50 800 600]);

subplot(2, 1, 1);
semilogx(vecteur_freq, 20*log10(spectre_positif + eps), ...
         'Color', [0.1 0.5 0.9], 'LineWidth', 0.8);
title('Spectre du Signal Original (échelle logarithmique)');
xlabel('Fréquence (Hz)');
ylabel('Magnitude (dB)');
xlim([10 freq_nyquist]);
grid on; grid minor;

subplot(2, 1, 2);
semilogx(vecteur_freq, 20*log10(spectre_filtre_final + eps), ...
         'Color', [0.9 0.1 0.1], 'LineWidth', 1.2);
title('Spectre du Signal Filtré (échelle logarithmique)');
xlabel('Fréquence (Hz)');
ylabel('Magnitude (dB)');
xlim([10 freq_nyquist]);
grid on; grid minor;

%% Partie 6: Analyse de la réponse du filtre global
fprintf('\n=== ÉTAPE 6: CARACTÉRISATION DU FILTRE ===\n');

% Calcul de la réponse en fréquence du filtre équivalent
nb_points_freq = 8192;
[reponse_complexe, freq_analyse] = freqz(numerateur_global, denominateur_global, ...
                                         nb_points_freq, freq_echantillonnage);

% Graphique 3: Caractéristiques du filtre
figure('Name', 'Réponse du filtre global', 'Position', [300 100 800 600]);

subplot(2, 1, 1);
semilogx(freq_analyse, 20*log10(abs(reponse_complexe)), ...
         'Color', [0.2 0.7 0.2], 'LineWidth', 1.5);
title('Gain du Filtre Global (Module)');
xlabel('Fréquence (Hz)');
ylabel('Gain (dB)');
grid on; grid minor;
xlim([1 freq_nyquist]);

subplot(2, 1, 2);
semilogx(freq_analyse, unwrap(angle(reponse_complexe)) * 180/pi, ...
         'Color', [0.7 0.2 0.7], 'LineWidth', 1.5);
title('Déphasage du Filtre Global');
xlabel('Fréquence (Hz)');
ylabel('Phase (degrés)');
grid on; grid minor;
xlim([1 freq_nyquist]);

%% Partie 7: Sauvegarde du résultat
fprintf('\n=== ÉTAPE 7: SAUVEGARDE ===\n');

% Nom du fichier de sortie
fichier_sortie = 'Mozart_Sans_Bruit.wav';

% Écriture du fichier audio filtré
audiowrite(fichier_sortie, signal_final, freq_echantillonnage, ...
           'BitsPerSample', bits_quantification);

fprintf('Signal nettoyé sauvegardé: %s\n', fichier_sortie);

%% Résumé des performances
fprintf('\n=== RÉSUMÉ DU TRAITEMENT ===\n');
puissance_originale = mean(signal_mozart.^2);
puissance_filtree = mean(signal_final.^2);
rapport_snr = 10*log10(puissance_filtree/puissance_originale);

fprintf('Nombre de filtres appliqués: %d\n', length(parasites_freq)-1);
fprintf('Ordre total du filtre: %d\n', length(numerateur_global)-1);
fprintf('Variation de puissance: %.2f dB\n', rapport_snr);
fprintf('Traitement terminé avec succès!\n');
