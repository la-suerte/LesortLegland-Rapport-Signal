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
% Date      : 14/05/2025


%% TP Traitement du Signal EEG - Filtrage par moyenne mobile
% Auteur: Étudiant 
% Date: Mai 2025
% Objectif: Appliquer un filtre à moyenne mobile sur des données EEG brutes

% Nettoyage de l'espace de travail
clear all; clc; close all;

%% Étape 1: Importation des données EEG
% Chargement du fichier contenant les mesures EEG
donnees_brutes = load('EEGDatabrut.mat');
signal_eeg = donnees_brutes.EEGDatabrut;

% Vérification de la taille du signal
fprintf('Nombre d''échantillons chargés: %d\n', length(signal_eeg));

%% Étape 2: Configuration des paramètres de traitement
freq_sampling = 100;                    % Fréquence d'échantillonnage (Hz)
taille_signal = length(signal_eeg);     % Nombre total d'échantillons
vecteur_temps = (0:taille_signal-1)/freq_sampling;  % Axe temporel en secondes

% Affichage des caractéristiques du signal
fprintf('Durée totale du signal: %.2f secondes\n', vecteur_temps(end));
fprintf('Résolution temporelle: %.4f secondes\n', 1/freq_sampling);

%% Étape 3: Implémentation du filtre à moyenne mobile (méthode manuelle)
% Création du vecteur de sortie initialisé à zéro
signal_filtre = zeros(size(signal_eeg));

% Traitement des premiers échantillons (conditions initiales)
% Pour n=1: y[1] = x[1] (pas d'échantillons précédents)
signal_filtre(1) = signal_eeg(1);

% Pour n=2: y[2] = (x[2] + x[1])/2
signal_filtre(2) = mean(signal_eeg(1:2));

% Pour n=3: y[3] = (x[3] + x[2] + x[1])/3
signal_filtre(3) = mean(signal_eeg(1:3));

% Pour n=4: y[4] = (x[4] + x[3] + x[2] + x[1])/4
signal_filtre(4) = mean(signal_eeg(1:4));

% Pour n=5: y[5] = (x[5] + x[4] + x[3] + x[2] + x[1])/5
signal_filtre(5) = mean(signal_eeg(1:5));

% Application de la formule du filtre à moyenne mobile pour n >= 6
% Équation: y[n] = (1/5) * somme(x[n-k]) pour k=0 à 4
ordre_filtre = 5;  % Nombre d'échantillons pour la moyenne
for indice = ordre_filtre+1:taille_signal
    % Calcul de la moyenne mobile sur 5 échantillons
    somme_echantillons = 0;
    for k = 0:ordre_filtre-1
        somme_echantillons = somme_echantillons + signal_eeg(indice-k);
    end
    signal_filtre(indice) = somme_echantillons / ordre_filtre;
end

fprintf('Filtrage terminé!\n');

%% Étape 4: Visualisation comparative des signaux
% Graphique 1: Vue d'ensemble sur toute la durée
figure('Name', 'Comparaison signal brut vs filtré');
plot(vecteur_temps, signal_eeg, 'Color', [0 0.4 0.8], 'LineWidth', 0.8);
hold on;
plot(vecteur_temps, signal_filtre, 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
xlabel('Temps (secondes)');
ylabel('Amplitude du signal EEG');
title('Comparaison entre signal EEG original et signal filtré');
legend('Signal EEG brut', 'Signal après filtrage', 'Location', 'best');
grid on;
set(gca, 'GridAlpha', 0.3);

% Graphique 2: Zoom sur les 2 premières secondes pour analyse détaillée
figure('Name', 'Analyse détaillée - 2 premières secondes');
plot(vecteur_temps, signal_eeg, 'Color', [0 0.6 0.9], 'LineWidth', 0.7);
hold on;
plot(vecteur_temps, signal_filtre, 'Color', [0.9 0.1 0.1], 'LineWidth', 2.2);
xlabel('Temps (secondes)');
ylabel('Amplitude');
title('Effet du filtrage sur les 2 premières secondes du signal EEG');
legend('EEG original', 'EEG filtré (moyenne mobile)', 'Location', 'northeast');
grid on;
xlim([0 2]);  % Limitation de l'affichage aux 2 premières secondes
set(gca, 'GridAlpha', 0.4);

%% Étape 5: Analyse spectrale du signal original
% Calcul de la transformée de Fourier discrète
nb_points = length(signal_eeg);
spectre_complexe = fft(signal_eeg);

% Calcul du spectre d'amplitude normalisé
spectre_amplitude = abs(spectre_complexe/nb_points);

% Extraction de la partie positive du spectre (0 à fe/2)
moitie_spectre = spectre_amplitude(1:floor(nb_points/2)+1);

% Correction pour les fréquences non-DC et non-Nyquist
moitie_spectre(2:end-1) = 2 * moitie_spectre(2:end-1);

% Création du vecteur de fréquences correspondant
step_freq = freq_sampling / nb_points;
axe_frequences = 0:step_freq:freq_sampling/2;

% Affichage du spectre fréquentiel
figure('Name', 'Analyse spectrale du signal EEG');
plot(axe_frequences, moitie_spectre, 'Color', [0.2 0.7 0.3], 'LineWidth', 1.5);
title('Spectre de fréquence du signal EEG brut');
xlabel('Fréquence (Hz)');
ylabel('Magnitude du spectre');
grid on;
xlim([0 50]);  % Affichage jusqu'à 50 Hz (fréquences physiologiques)
set(gca, 'GridAlpha', 0.3);

% Affichage de statistiques sur le signal
fprintf('\n=== STATISTIQUES DU TRAITEMENT ===\n');
fprintf('Valeur moyenne du signal original: %.4f\n', mean(signal_eeg));
fprintf('Valeur moyenne du signal filtré: %.4f\n', mean(signal_filtre));
fprintf('Écart-type du signal original: %.4f\n', std(signal_eeg));
fprintf('Écart-type du signal filtré: %.4f\n', std(signal_filtre));
fprintf('Réduction du bruit: %.2f%%\n', (1-std(signal_filtre)/std(signal_eeg))*100);
