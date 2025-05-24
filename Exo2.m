%%EX2
%% Filtrage d'un bruit aigu dans un signal audio
% Ce script analyse un fichier audio pour identifier et
% supprimer des bruits aigus indésirables
% Nettoyage de l'environnement
clc; clear; close all;
% 1. Chargement des données
[melodie_signal, freq_echant] = audioread('Mozart_Bruit.wav');
info_audio = audioinfo('Mozart_Bruit.wav');
% Paramètres
duree_totale = length(melodie_signal) / freq_echant;
nb_bits = info_audio.BitsPerSample;
nb_echant = length(melodie_signal);
axe_temps = (0:nb_echant-1) / freq_echant;
% Affichage des caractéristiques
fprintf('Période d''échantillonnage: %.6f secondes\n', 1/freq_echant);
fprintf('Fréquence d''échantillonnage: %d Hz\n', freq_echant);
fprintf('Nombre de bits pour le codage: %d bits\n', nb_bits);
fprintf('Durée totale du signal: %.2f secondes\n', duree_totale);
% Écouter le son (décommenter pour écouter)
% sound(melodie_signal, freq_echant);
% Calcul de la transformée de Fourier
transformee_fourier = fft(melodie_signal);
spectre_ampli = abs(transformee_fourier/nb_echant);
spectre_ampli_mono = spectre_ampli(1:floor(nb_echant/2)+1);
spectre_ampli_mono(2:end-1) = 2*spectre_ampli_mono(2:end-1);
% Vecteur de fréquence pour l'affichage
resolution_freq = freq_echant/nb_echant;
axe_freq = 0:resolution_freq:freq_echant/2;
% Conversion en dB et limiter pour éviter log(0)
epsilon = 1e-10;
spectre_db = 20*log10(spectre_ampli_mono + epsilon);
% Figure 1: Représentation temporelle et fréquentielle
figure;
subplot(2,1,1);
plot(axe_temps, melodie_signal, 'b', 'LineWidth', 0.5);
xlabel('Temps (s)');
ylabel('Amplitude');
title('Représentation temporelle du signal');
grid on;
% Représentation fréquentielle
subplot(2,1,2);
plot(axe_freq, spectre_db, 'r', 'LineWidth', 1);
xlabel('Fréquence (Hz)');
ylabel('Amplitude (dB)');
title('Spectre d''amplitude en dB');
grid on;
xlim([0 freq_echant/2]);
% Détermination des fréquences indésirables
[valeurs_pics, positions_pics] = findpeaks(spectre_db, 'MinPeakHeight', max(spectre_db)-20, 'MinPeakDistance', 500);
freq_pics = axe_freq(positions_pics);
% Figure 2: Identification des fréquences indésirables
figure;
plot(axe_freq, spectre_db, 'b', 'LineWidth', 1);
hold on;
plot(freq_pics, valeurs_pics, 'ro', 'MarkerSize', 10);
title('Spectre avec identification des fréquences indésirables');
xlabel('Fréquence (Hz)');
ylabel('Amplitude (dB)');
grid on;
xlim([0 freq_echant/2]);
legend('Spectre', 'Pics détectés');
% Affichage des valeurs des fréquences indésirables détectées
fprintf('Fréquences indésirables détectées :\n');
for indice = 1:length(freq_pics)
fprintf('Pic %d: %.2f Hz\n', indice, freq_pics(indice));
end
% Spectre avec fftshift pour visualisation [-fe/2, fe/2]
spectre_complet = abs(fftshift(transformee_fourier))/nb_echant;
spectre_complet_db = 20*log10(spectre_complet + epsilon);
% Axe des fréquences centré
axe_freq_centre = linspace(-freq_echant/2, freq_echant/2, nb_echant);
% Figure 3: Représentation du spectre centré
figure;
plot(axe_freq_centre, spectre_complet_db, 'g', 'LineWidth', 1);
title('Spectre d''amplitude en dB sur [-f_e/2, f_e/2]');
xlabel('Fréquence (Hz)');
ylabel('Amplitude (dB)');
grid on;
xlim([-freq_echant/2 freq_echant/2]);
% Détermination des fréquences indésirables dans l'intervalle [0, fe/2]
[valeurs_pics_centre, positions_pics_centre] = findpeaks(spectre_complet_db(floor(nb_echant/2):end), 'MinPeakHeight', max(spectre_complet_db)-30, 'MinPeakDistance', 500);
freq_pics_centre = axe_freq_centre(positions_pics_centre + floor(nb_echant/2) - 1);
% Filtrer pour n'afficher que les fréquences entre 0 et fe/2
freq_pics_positives = freq_pics_centre(freq_pics_centre > 0 & freq_pics_centre < freq_echant/2);
fprintf('\nFréquences indésirables entre 0 et fe/2 :\n');
for indice = 1:length(freq_pics_positives)
fprintf('Pic %d: %.2f Hz\n', indice, freq_pics_positives(indice));
end
