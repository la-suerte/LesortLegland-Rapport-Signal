%%Exercice 1
%    ______  _____ _____ _  __      __  
%   |  ____|/ ____|_   _| | \ \    / /  
%   | |__  | (___   | | | |  \ \  / /  
%   |  __|  \___ \  | | | |   \ \/ /    
%   | |____ ____) |_| | | |____\  /  
%   |______|_____/|_____|______|\/      
%                                      
%   E N G I N E E R I N G   S C H O O L
%   D E   V I N C I   P A R I S        
%
% Copyright (c) 2025. All rights reserved.
%
% Subject   : Signal Processing
% Project   : TP Matlab
% Author(s) : François-Louis Legland & Antoine G. Lesort
% Team      : J
% Date      : 24/05/2025
%% TP5 - Démodulation FSK (Frequency Shift Keying)
% Objectif: Décoder des messages binaires transmis par modulation de fréquence
% Méthode: Utilisation de filtres passe-bande sélectifs pour séparer f1 et f2
% Auteur: Étudiant - Communications Numériques

% Initialisation de l'environnement
clear variables; clc; close all;

fprintf('=== DÉMODULATION DE SIGNAUX FSK ===\n\n');

%% Configuration des paramètres généraux
freq_echantillonnage = 100e3;      % Fréquence d'échantillonnage (100 kHz)
duree_symbole = 1e-3;              % Durée d'un symbole Ts = 1 ms
periode_ech = 1/freq_echantillonnage;  % Période d'échantillonnage

fprintf('Paramètres du système FSK:\n');
fprintf('- Fréquence d''échantillonnage: %.0f kHz\n', freq_echantillonnage/1000);
fprintf('- Durée d''un symbole: %.1f ms\n', duree_symbole*1000);
fprintf('- Résolution temporelle: %.2f µs\n', periode_ech*1e6);

%% Fonction pour créer un filtre passe-bande de 2ème ordre
function [coef_num, coef_den] = creer_filtre_pb(f_centrale, r_val, fe)
    % Calcul des coefficients pour H(z) = 1/(1 - 2r*cos(2πf0*Te)*z^(-1) + r^2*z^(-2))
    Te = 1/fe;
    freq_norm = f_centrale * Te;  % Fréquence normalisée
    
    % Coefficients du dénominateur
    coef_z1 = 2 * r_val * cos(2*pi*freq_norm);
    coef_z2 = r_val^2;
    
    coef_num = 1;                           % Numérateur constant
    coef_den = [1, -coef_z1, coef_z2];     % Dénominateur
end

%% Fonction pour analyser et décoder un signal FSK
function decoder_signal_fsk(nom_fichier, f1_hz, f2_hz, fe, Ts)
    fprintf('\n=== ANALYSE DU SIGNAL: %s ===\n', nom_fichier);
    fprintf('Fréquences FSK: f1 = %.0f kHz (bit 1), f2 = %.0f kHz (bit 0)\n', ...
            f1_hz/1000, f2_hz/1000);
    
    %% Étape 1: Chargement et affichage du signal
    donnees_signal = load(nom_fichier);
    noms_champs = fieldnames(donnees_signal);
    signal_fsk = donnees_signal.(noms_champs{1});  % Premier champ du fichier
    
    nb_echantillons = length(signal_fsk);
    duree_totale = nb_echantillons / fe;
    vecteur_temps = (0:nb_echantillons-1) / fe;
    
    fprintf('\nCaractéristiques du signal chargé:\n');
    fprintf('- Nombre d''échantillons: %d\n', nb_echantillons);
    fprintf('- Durée totale: %.2f ms\n', duree_totale*1000);
    fprintf('- Nombre de symboles estimé: %.0f\n', duree_totale/Ts);
    
    % Visualisation temporelle
    figure('Name', sprintf('Signal FSK - %s', nom_fichier), ...
           'Position', [50 + rand*100, 50 + rand*100, 1200, 800]);
    
    subplot(3, 2, 1);
    plot(vecteur_temps*1000, signal_fsk, 'b-', 'LineWidth', 1);
    title(sprintf('Signal FSK Reçu - %s', nom_fichier));
    xlabel('Temps (ms)');
    ylabel('Amplitude');
    grid on; grid minor;
    
    %% Étape 2: Analyse spectrale avec FFT
    spectre_fft = fft(signal_fsk);
    magnitude_spectre = abs(spectre_fft) / nb_echantillons;
    
    % Construction du spectre unilatéral
    spectre_positif = magnitude_spectre(1:floor(nb_echantillons/2)+1);
    spectre_positif(2:end-1) = 2 * spectre_positif(2:end-1);
    
    % Vecteur des fréquences
    resolution_freq = fe / nb_echantillons;
    frequences = (0:floor(nb_echantillons/2)) * resolution_freq;
    
    subplot(3, 2, 2);
    plot(frequences/1000, 20*log10(spectre_positif + eps), 'r-', 'LineWidth', 1.5);
    title('Spectre Fréquentiel du Signal FSK');
    xlabel('Fréquence (kHz)');
    ylabel('Magnitude (dB)');
    grid on; grid minor;
    xlim([0 50]);  % Affichage jusqu'à 50 kHz
    
    % Identification des pics spectraux
    [pics_valeurs, pics_positions] = findpeaks(spectre_positif, ...
        'MinPeakHeight', max(spectre_positif)*0.1, 'MinPeakDistance', 100);
    frequences_pics = frequences(pics_positions);
    
    fprintf('\nPics spectraux détectés:\n');
    for i = 1:length(frequences_pics)
        fprintf('  - Pic %d: %.1f kHz\n', i, frequences_pics(i)/1000);
    end
    
    %% Étape 3: Conception des filtres passe-bande
    facteur_r = 0.9;  % Facteur de sélectivité
    
    % Filtre pour f1 (détection des '1')
    [num_f1, den_f1] = creer_filtre_pb(f1_hz, facteur_r, fe);
    
    % Filtre pour f2 (détection des '0')
    [num_f2, den_f2] = creer_filtre_pb(f2_hz, facteur_r, fe);
    
    % Calcul des réponses fréquentielles
    nb_points_freq = 2048;
    [reponse_f1, freq_analyse] = freqz(num_f1, den_f1, nb_points_freq, fe);
    [reponse_f2, ~] = freqz(num_f2, den_f2, nb_points_freq, fe);
    
    subplot(3, 2, 3);
    plot(freq_analyse/1000, 20*log10(abs(reponse_f1)), 'g-', 'LineWidth', 2);
    hold on;
    plot(freq_analyse/1000, 20*log10(abs(reponse_f2)), 'm-', 'LineWidth', 2);
    title('Réponses Fréquentielles des Filtres');
    xlabel('Fréquence (kHz)');
    ylabel('Gain (dB)');
    legend(sprintf('Filtre f1 (%.0f kHz)', f1_hz/1000), ...
           sprintf('Filtre f2 (%.0f kHz)', f2_hz/1000), 'Location', 'best');
    grid on; grid minor;
    xlim([0 50]);
    
    %% Étape 4: Filtrage du signal FSK
    % Application des filtres
    signal_filtre_f1 = filter(num_f1, den_f1, signal_fsk);
    signal_filtre_f2 = filter(num_f2, den_f2, signal_fsk);
    
    % Calcul des enveloppes (détection d'amplitude)
    enveloppe_f1 = abs(hilbert(signal_filtre_f1));
    enveloppe_f2 = abs(hilbert(signal_filtre_f2));
    
    subplot(3, 2, 4);
    plot(vecteur_temps*1000, enveloppe_f1, 'g-', 'LineWidth', 1.5);
    hold on;
    plot(vecteur_temps*1000, enveloppe_f2, 'm-', 'LineWidth', 1.5);
    title('Enveloppes des Signaux Filtrés');
    xlabel('Temps (ms)');
    ylabel('Amplitude');
    legend('Enveloppe f1 (bit 1)', 'Enveloppe f2 (bit 0)', 'Location', 'best');
    grid on; grid minor;
    
    %% Étape 5: Décodage binaire
    % Échantillonnage au milieu de chaque symbole
    nb_symboles = floor(duree_totale / Ts);
    echant_par_symbole = round(Ts * fe);
    
    sequence_binaire = zeros(1, nb_symboles);
    
    for k = 1:nb_symboles
        % Position du milieu du k-ième symbole
        pos_echantillon = round((k-0.5) * echant_par_symbole);
        
        if pos_echantillon > 0 && pos_echantillon <= nb_echantillons
            % Décision basée sur les enveloppes
            energie_f1 = enveloppe_f1(pos_echantillon);
            energie_f2 = enveloppe_f2(pos_echantillon);
            
            if energie_f1 > energie_f2
                sequence_binaire(k) = 1;
            else
                sequence_binaire(k) = 0;
            end
        end
    end
    
    % Visualisation de la décision
    temps_symboles = ((1:nb_symboles) - 0.5) * Ts * 1000;  % en ms
    
    subplot(3, 2, 5);
    stem(temps_symboles, sequence_binaire, 'ko-', 'LineWidth', 2, 'MarkerSize', 8);
    title('Séquence Binaire Décodée');
    xlabel('Temps (ms)');
    ylabel('Bit Décodé');
    ylim([-0.5 1.5]);
    grid on; grid minor;
    
    %% Étape 6: Conversion en décimal
    % Affichage de la séquence binaire
    chaine_binaire = sprintf('%d', sequence_binaire);
    fprintf('\nRésultats du décodage:\n');
    fprintf('- Séquence binaire: %s\n', chaine_binaire);
    fprintf('- Longueur: %d bits\n', length(sequence_binaire));
    
    % Conversion en décimal (MSB en premier)
    valeur_decimale = 0;
    for i = 1:length(sequence_binaire)
        valeur_decimale = valeur_decimale + sequence_binaire(i) * 2^(length(sequence_binaire)-i);
    end
    
    fprintf('- Valeur décimale: %d\n', valeur_decimale);
    
    % Conversion par groupes de 8 bits (octets)
    if length(sequence_binaire) >= 8
        nb_octets = floor(length(sequence_binaire) / 8);
        fprintf('- Conversion par octets:\n');
        
        for octet = 1:nb_octets
            debut = (octet-1)*8 + 1;
            fin = octet*8;
            bits_octet = sequence_binaire(debut:fin);
            
            valeur_octet = 0;
            for bit = 1:8
                valeur_octet = valeur_octet + bits_octet(bit) * 2^(8-bit);
            end
            
            chaine_octet = sprintf('%d', bits_octet);
            fprintf('    Octet %d: %s = %d\n', octet, chaine_octet, valeur_octet);
        end
    end
    
    % Histogramme des décisions
    subplot(3, 2, 6);
    histogram(sequence_binaire, [-0.5, 0.5, 1.5], 'FaceColor', [0.7 0.7 0.9]);
    title('Distribution des Bits Décodés');
    xlabel('Valeur du Bit');
    ylabel('Nombre d''Occurrences');
    xticks([0 1]);
    xticklabels({'0', '1'});
    grid on;
    
    fprintf('- Nombre de ''0'': %d\n', sum(sequence_binaire == 0));
    fprintf('- Nombre de ''1'': %d\n', sum(sequence_binaire == 1));
end

%% ANALYSE DU PREMIER SIGNAL
decoder_signal_fsk('signal.mat', 40e3, 20e3, freq_echantillonnage, duree_symbole);

%% ANALYSE DU SECOND SIGNAL
decoder_signal_fsk('signal2.mat', 35e3, 40e3, freq_echantillonnage, duree_symbole);

%% Synthèse et comparaison
fprintf('\n\n=== SYNTHÈSE DE L''ANALYSE FSK ===\n');
fprintf('Principe de démodulation FSK:\n');
fprintf('1. Analyse spectrale pour identifier f1 et f2\n');
fprintf('2. Filtrage sélectif avec filtres passe-bande\n');
fprintf('3. Détection d''enveloppe pour mesurer l''énergie\n');
fprintf('4. Échantillonnage au centre de chaque symbole\n');
fprintf('5. Décision binaire basée sur les énergies relatives\n\n');

fprintf('Paramètres des filtres utilisés:\n');
fprintf('- Ordre: 2 (filtre récursif)\n');
fprintf('- Facteur r: 0.9 (haute sélectivité)\n');
fprintf('- Type: Passe-bande centré sur f1 et f2\n\n');

fprintf('Avantages de cette méthode:\n');
fprintf('- Robuste au bruit\n');
fprintf('- Implémentation simple\n');
fprintf('- Faible complexité calculatoire\n');

fprintf('\nAnalyse FSK terminée!\n');
