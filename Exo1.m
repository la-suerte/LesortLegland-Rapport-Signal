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
%% Chargement des données
data = load('EEGDatabrut.mat');
x = data.EEGDatabrut;
%% Paramètres
fe = 100; % Fréquence d'échantillonnage en Hz
N = length(x); % Longueur du signal
t = (0:N-1)/fe; % Vecteur temps
% Méthode 1: Implémentation manuelle de l'équation aux différences pour le filtre à moyenne mobile
y_manual = zeros(1, N);
%% Conditions initiales (premiers échantillons)
y_manual(1) = x(1);
y_manual(2) = (x(2) + x(1))/2;
y_manual(3) = (x(3) + x(2) + x(1))/3;
y_manual(4) = (x(4) + x(3) + x(2) + x(1))/4;
y_manual(5) = (x(5) + x(4) + x(3) + x(2) + x(1))/5;
% Boucle principale pour calculer la sortie avec l'équation demandée
for n = 5:N
    y_manual(n) = (1/5) * (x(n) + x(n-1) + x(n-2) + x(n-3) + x(n-4));
end
%% Figure 1: Signal brut et filtré sur toute la durée
figure;
plot(t, x, 'b', 'LineWidth', 0.5);
hold on;
plot(t, y_manual, 'r', 'LineWidth', 1.5);
xlabel('Temps (s)');
ylabel('Amplitude');
title('Signal EEG brut et filtré (moyenne mobile)');
legend('Signal brut', 'Signal filtré');
grid on;
%% Figure 2: Zoom sur les deux premières secondes pour le comptage des pics
figure;
plot(t, x, 'b', 'LineWidth', 0.5);
hold on;
plot(t, y_manual, 'r', 'LineWidth', 1.5);
xlabel('Temps (s)');
ylabel('Amplitude');
title('Signal EEG brut et filtré (2 premières secondes)');
legend('Signal brut', 'Signal filtré');
grid on;
xlim([0 2]);
%% Calcul et affichage du spectre du signal brut
N = length(x);
Y = fft(x);
P2 = abs(Y/N);
P1 = P2(1:floor(N/2)+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fe*(0:(N/2))/N;
figure;
plot(f, P1);
title('Spectre du signal EEG brut');
xlabel('Fréquence (Hz)');
ylabel('Amplitude');
grid on;
xlim([0 50]); % Affichage jusqu'à 50 Hz
