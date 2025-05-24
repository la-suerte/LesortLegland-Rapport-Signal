%%Exercice 4
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
%% TP4 - Analyse d'un Filtre Récursif de Second Ordre
% Fonction de transfert: H(z) = 1 / (1 - 2r*cos(2πf0*Te)*z^(-1) + r^2*z^(-2))
% Objectif: Étudier l'influence des paramètres r et f0 sur les caractéristiques du filtre
% Auteur: Étudiant - Traitement Numérique du Signal

% Nettoyage de l'espace de travail
clear all; clc; close all;

%% Configuration générale des paramètres
fprintf('=== ANALYSE DU FILTRE RÉCURSIF D''ORDRE 2 ===\n\n');

% Paramètres de base
freq_sampling = 1;        % Fréquence d'échantillonnage (Hz)
periode_ech = 1/freq_sampling;  % Période d'échantillonnage (s)
nb_echantillons = 512;    % Nombre de points pour l'analyse
duree_analyse = 50;       % Durée pour les réponses temporelles

fprintf('Paramètres d''analyse:\n');
fprintf('- Fréquence d''échantillonnage: %.1f Hz\n', freq_sampling);
fprintf('- Période d''échantillonnage: %.1f s\n', periode_ech);
fprintf('- Points d''analyse fréquentielle: %d\n', nb_echantillons);

%% Fonction pour analyser un filtre donné
function analyser_filtre(r_val, f0_val, fe, titre_config)
    % Calcul des coefficients du dénominateur
    Te = 1/fe;
    coef_z1 = 2 * r_val * cos(2*pi*f0_val*Te);  % Coefficient de z^(-1)
    coef_z2 = r_val^2;                          % Coefficient de z^(-2)
    
    % Définition des coefficients du filtre
    numerateur = 1;                             % Numérateur constant
    denominateur = [1, -coef_z1, coef_z2];     % 1 - coef_z1*z^(-1) + coef_z2*z^(-2)
    
    fprintf('\n--- %s ---\n', titre_config);
    fprintf('Paramètres: r = %.2f, f0 = %.3f Hz\n', r_val, f0_val);
    fprintf('Coefficients du dénominateur: [1, %.4f, %.4f]\n', -coef_z1, coef_z2);
    
    % Vérification de la stabilité
    poles_filtre = roots(denominateur);
    module_poles = abs(poles_filtre);
    stabilite = all(module_poles < 1);
    
    fprintf('Pôles du filtre:\n');
    for i = 1:length(poles_filtre)
        fprintf('  Pôle %d: %.4f%+.4fi (module = %.4f)\n', i, ...
                real(poles_filtre(i)), imag(poles_filtre(i)), module_poles(i));
    end
    fprintf('Stabilité: %s\n', char("STABLE" * stabilite + "INSTABLE" * ~stabilite));
    
    % Calcul de la réponse en fréquence
    [reponse_freq, freq_vect] = freqz(numerateur, denominateur, nb_echantillons, fe);
    magnitude_db = 20*log10(abs(reponse_freq));
    phase_deg = angle(reponse_freq) * 180/pi;
    
    % Réponse impulsionnelle
    duree_imp = 40;
    reponse_impuls = impz(numerateur, denominateur, duree_imp);
    temps_imp = 0:duree_imp-1;
    
    % Réponse indicielle (échelon unité)
    signal_echelon = ones(duree_imp, 1);
    reponse_echelon = filter(numerateur, denominateur, signal_echelon);
    
    % Création des graphiques
    figure('Name', sprintf('Analyse - r=%.2f, f0=%.3f', r_val, f0_val), ...
           'Position', [100 + rand*200, 100 + rand*200, 1000, 700]);
    
    % Sous-graphique 1: Réponse en fréquence (magnitude)
    subplot(2, 3, 1);
    plot(freq_vect, magnitude_db, 'b-', 'LineWidth', 2);
    title(sprintf('Réponse en Fréquence - Magnitude\nr=%.2f, f_0=%.3f Hz', r_val, f0_val));
    xlabel('Fréquence (Hz)');
    ylabel('Gain (dB)');
    grid on; grid minor;
    ylim([-60 20]);
    
    % Sous-graphique 2: Réponse en fréquence (phase)
    subplot(2, 3, 2);
    plot(freq_vect, phase_deg, 'r-', 'LineWidth', 2);
    title('Phase');
    xlabel('Fréquence (Hz)');
    ylabel('Phase (degrés)');
    grid on; grid minor;
    
    % Sous-graphique 3: Plan des pôles et zéros
    subplot(2, 3, 3);
    % Cercle unité pour référence
    theta = 0:0.01:2*pi;
    plot(cos(theta), sin(theta), 'k--', 'LineWidth', 1);
    hold on;
    % Affichage des pôles
    plot(real(poles_filtre), imag(poles_filtre), 'rx', 'MarkerSize', 10, 'LineWidth', 3);
    % Zéro à l'origine (implicite)
    plot(0, 0, 'bo', 'MarkerSize', 8, 'LineWidth', 2);
    title('Pôles et Zéros');
    xlabel('Partie Réelle');
    ylabel('Partie Imaginaire');
    grid on; axis equal;
    xlim([-1.5 1.5]); ylim([-1.5 1.5]);
    legend('Cercle unité', 'Pôles', 'Zéros', 'Location', 'best');
    
    % Sous-graphique 4: Réponse impulsionnelle
    subplot(2, 3, 4);
    stem(temps_imp, reponse_impuls, 'g-', 'LineWidth', 1.5);
    title('Réponse Impulsionnelle');
    xlabel('Échantillons (n)');
    ylabel('Amplitude h[n]');
    grid on; grid minor;
    
    % Sous-graphique 5: Réponse indicielle
    subplot(2, 3, 5);
    stem(temps_imp, reponse_echelon, 'm-', 'LineWidth', 1.5);
    title('Réponse Indicielle');
    xlabel('Échantillons (n)');
    ylabel('Amplitude s[n]');
    grid on; grid minor;
    
    % Sous-graphique 6: Diagramme de Bode
    subplot(2, 3, 6);
    semilogx(freq_vect, magnitude_db, 'k-', 'LineWidth', 2);
    title('Diagramme de Bode');
    xlabel('Fréquence (Hz)');
    ylabel('Gain (dB)');
    grid on; grid minor;
    xlim([0.001 0.5]);
    
    % Identification du type de filtre
    [gain_max, idx_max] = max(magnitude_db);
    freq_resonance = freq_vect(idx_max);
    
    % Analyse du comportement
    gain_dc = magnitude_db(1);  % Gain à f=0
    gain_nyquist = magnitude_db(end);  % Gain à f=fe/2
    
    fprintf('Caractéristiques du filtre:\n');
    fprintf('  - Gain DC (f=0): %.2f dB\n', gain_dc);
    fprintf('  - Gain max: %.2f dB à f=%.4f Hz\n', gain_max, freq_resonance);
    fprintf('  - Gain Nyquist (f=fe/2): %.2f dB\n', gain_nyquist);
    
    % Classification du type de filtre
    if abs(freq_resonance) < 0.01  % Proche de f=0
        type_filtre = 'Passe-bas';
    elseif abs(freq_resonance - 0.5) < 0.01  % Proche de fe/2
        type_filtre = 'Passe-haut';
    elseif freq_resonance > 0.01 && freq_resonance < 0.49
        type_filtre = 'Passe-bande';
    else
        type_filtre = 'Indéterminé';
    end
    
    fprintf('  - Type de filtre identifié: %s\n', type_filtre);
end

%% PARTIE 1: Étude de la stabilité avec f0 = 0.1 Hz
fprintf('\n\n=== PARTIE 1: ÉTUDE DE LA STABILITÉ (f0 = 0.1 Hz) ===\n');
f0_test = 0.1;
valeurs_r = 0:0.5:2;

for r_actuel = valeurs_r
    analyser_filtre(r_actuel, f0_test, freq_sampling, ...
                   sprintf('Test stabilité - r=%.1f', r_actuel));
    pause(1); % Pause pour voir les résultats
end

%% PARTIE 2: Filtre avec f0 = 0 Hz (Passe-bas)
fprintf('\n\n=== PARTIE 2: COMPORTEMENT PASSE-BAS (f0 = 0 Hz) ===\n');
f0_test = 0;
valeurs_r = [0.1, 0.3, 0.5, 0.7, 0.9];

for r_actuel = valeurs_r
    analyser_filtre(r_actuel, f0_test, freq_sampling, ...
                   sprintf('Passe-bas - r=%.1f', r_actuel));
    pause(1);
end

%% PARTIE 3: Filtre avec f0 = 0.5 Hz (Passe-haut)
fprintf('\n\n=== PARTIE 3: COMPORTEMENT PASSE-HAUT (f0 = 0.5 Hz) ===\n');
f0_test = 0.5;
valeurs_r = [0.1, 0.3, 0.5, 0.7, 0.9];

for r_actuel = valeurs_r
    analyser_filtre(r_actuel, f0_test, freq_sampling, ...
                   sprintf('Passe-haut - r=%.1f', r_actuel));
    pause(1);
end

%% PARTIE 4: Filtre avec f0 = 0.25 Hz (Passe-bande)
fprintf('\n\n=== PARTIE 4: COMPORTEMENT PASSE-BANDE (f0 = 0.25 Hz) ===\n');
f0_test = 0.25;
valeurs_r = [0.1, 0.3, 0.5, 0.7, 0.9];

for r_actuel = valeurs_r
    analyser_filtre(r_actuel, f0_test, freq_sampling, ...
                   sprintf('Passe-bande - r=%.1f', r_actuel));
    pause(1);
end

%% Synthèse théorique
fprintf('\n\n=== SYNTHÈSE THÉORIQUE ===\n');
fprintf('Expression des pôles:\n');
fprintf('  p1,2 = r * exp(±j*2π*f0*Te)\n');
fprintf('  En coordonnées cartésiennes: r*cos(2πf0Te) ± j*r*sin(2πf0Te)\n\n');

fprintf('Condition de stabilité:\n');
fprintf('  |p1| < 1 et |p2| < 1\n');
fprintf('  Soit: r < 1 (condition nécessaire et suffisante)\n\n');

fprintf('Influence des paramètres:\n');
fprintf('  - r contrôle la proximité des pôles au cercle unité\n');
fprintf('  - f0 détermine la fréquence de résonance du filtre\n');
fprintf('  - Plus r est proche de 1, plus le filtre est sélectif\n\n');

fprintf('Types de filtres obtenus:\n');
fprintf('  - f0 = 0      → Filtre passe-bas\n');
fprintf('  - f0 = fe/2   → Filtre passe-haut\n');
fprintf('  - 0 < f0 < fe/2 → Filtre passe-bande centré sur f0\n');

fprintf('\nAnalyse terminée!\n');
