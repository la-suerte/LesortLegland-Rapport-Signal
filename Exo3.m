%%EX3

% Filtrage numérique en temps réel avec filtre IIR causal
% Ce programme permet de filtrer un signal audio pour supprimer des bruits indésirables
% Nettoyage de l'environnement
clc; clear; close all;
% 1. Chargement des données
cheminFichier = 'Mozart_Bruit.wav';
[audioSource, tauxEchan] = audioread(cheminFichier);
% Conversion en mono si le signal est stéréo
if size(audioSource, 2) > 1
audioSource = mean(audioSource, 2);
end
% Paramètres
info_audio = audioinfo(cheminFichier);
dureeSignal = length(audioSource) / tauxEchan;
periodEchan = 1 / tauxEchan;
nbBits = info_audio.BitsPerSample;
longSignal = length(audioSource);
axeTemps = (0:longSignal-1) * periodEchan;
% Affichage des caractéristiques
fprintf('Période d''échantillonnage: %.6f secondes\n', periodEchan);
fprintf('Fréquence d''échantillonnage: %d Hz\n', tauxEchan);
fprintf('Nombre de bits pour le codage: %d bits\n', nbBits);
fprintf('Durée totale du signal: %.2f secondes\n', dureeSignal);
% Calcul de la transformée de Fourier du signal d'origine
axeFreq = (0:longSignal/2) * tauxEchan / longSignal;
spectreFourier = abs(fft(audioSource)) / longSignal;
spectreFourier = spectreFourier(1:longSignal/2+1);
spectreFourier(2:end-1) = 2 * spectreFourier(2:end-1);
% Création des filtres Butterworth pour chaque fréquence indésirable
freqNyquist = tauxEchan / 2;
audioFiltre = audioSource;
% Liste des fréquences indésirables à filtrer
freqIndesirables = [0.03, 117.01, 160.93, 333.19, 1397.01, 5658.01];
coefsTotNum = 1;
coefsTotDen = 1;
% Application des filtres en cascade
for i = 1:length(freqIndesirables)
freqActuelle = freqIndesirables(i);
% Ignorer les très basses fréquences
if freqActuelle < 1
continue;
end
% Calcul des limites de la bande à couper (±25 Hz)
freqBasse = max(0, (freqActuelle - 25) / freqNyquist);
freqHaute = min(1, (freqActuelle + 25) / freqNyquist);
% Création du filtre coupe-bande
[coefsNum, coefsDen] = butter(2, [freqBasse, freqHaute], 'stop');
% Filtrage du signal
audioFiltre = filter(coefsNum, coefsDen, audioFiltre);
% Mise à jour des coefficients du filtre total
coefsTotNum = conv(coefsTotNum, coefsNum);
coefsTotDen = conv(coefsTotDen, coefsDen);
fprintf('Filtre %d: Fréquence centrale = %.2f Hz, Bande coupée = [%.2f, %.2f] Hz\n', ...
i, freqActuelle, freqActuelle-25, freqActuelle+25);
end
% Normalisation du signal filtré
facteurMax = max(abs(audioFiltre));
audioNormalise = audioFiltre / facteurMax;
% Calcul du spectre du signal filtré
spectreFiltreFFT = abs(fft(audioNormalise)) / longSignal;
spectreFiltreFFT = spectreFiltreFFT(1:longSignal/2+1);
spectreFiltreFFT(2:end-1) = 2 * spectreFiltreFFT(2:end-1);
% Figure 1: Représentation temporelle
figure;
subplot(2, 1, 1);
plot(axeTemps, audioSource, 'b', 'LineWidth', 0.5);
title('Signal Audio Original');
xlabel('Temps (s)');
ylabel('Amplitude');
grid on;
subplot(2, 1, 2);
plot(axeTemps, audioNormalise, 'r', 'LineWidth', 1);
title('Signal Audio Filtré');
xlabel('Temps (s)');
ylabel('Amplitude');
grid on;
% Figure 2: Représentation fréquentielle
figure;
subplot(2, 1, 1);
semilogx(axeFreq, 20*log10(spectreFourier + eps), 'b', 'LineWidth', 0.5);
title('Spectre du Signal Original');
xlabel('Fréquence (Hz)');
ylabel('Amplitude (dB)');
xlim([10 freqNyquist]);
grid on;
subplot(2, 1, 2);
semilogx(axeFreq, 20*log10(spectreFiltreFFT + eps), 'r', 'LineWidth', 1);
title('Spectre du Signal Filtré');
xlabel('Fréquence (Hz)');
ylabel('Amplitude (dB)');
xlim([10 freqNyquist]);
grid on;
% Sauvegarde du signal filtré
cheminSortie = 'Mozart_Filtre.wav';
audiowrite(cheminSortie, audioNormalise, tauxEchan, 'BitsPerSample', nbBits);
%% Figure 3: Réponse en fréquence du filtre équivalent
nbPoints = 8192;
[reponseFreq, vectFreq] = freqz(coefsTotNum, coefsTotDen, nbPoints, tauxEchan);
figure;
subplot(2, 1, 1);
semilogx(vectFreq, 20*log10(abs(reponseFreq)), 'g', 'LineWidth', 1);
title('Module de la Réponse en Fréquence');
xlabel('Fréquence (Hz)');
ylabel('Gain (dB)');
grid on;
xlim([10 freqNyquist]);
subplot(2, 1, 2);
semilogx(vectFreq, unwrap(angle(reponseFreq)) * 180/pi, 'g', 'LineWidth', 1);
title('Phase de la Réponse en Fréquence');
xlabel('Fréquence (Hz)');
ylabel('Phase (degrés)');
grid on;
xlim([10 freqNyquist]);
% Affichage des informations sur le filtre total
fprintf('\nTraitement terminé. Signal filtré enregistré sous %s\n', cheminSortie);
