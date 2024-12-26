clear all
close all
clc

data = readmatrix("DecDec2.txt", "Delimiter", ",");
% NON VA

%%
% NON VA
% Apri il file
fid = fopen('Dec24_14_15.txt', 'r');

% Leggi i dati come numeri separati da virgole
dataArray = textscan(fid, repmat('%f', 1, 12), 'Delimiter', ',');

% Chiudi il file
fclose(fid);

% Combina i dati in una matrice
data2 = cell2mat(dataArray);

% Visual

%% FUNZIONA
% Apri il file
fid = fopen('DecDec2.txt', 'r');

% Leggi i dati come stringhe di caratteri
textData = textscan(fid, '%s', 'Delimiter', ';');

% Chiudi il file
fclose(fid);

% Rimuovi eventuali spazi extra
cleanedData = strrep(textData{1}, ' ', '');

% Inizializza una matrice vuota per i dati
data = [];

% Converti le stringhe di caratteri in numeri
for i = 1:length(cleanedData)
    numArray = str2double(strsplit(cleanedData{i}, ','));
    data = [data; numArray];
end

% Visualizza la matrice dei dati
disp(data);

%%
% NON VA
% Apri il file
fid = fopen('DecDec2.txt', 'r');

% Leggi tutto il contenuto del file come testo
fileContent = fread(fid, '*char')';

% Chiudi il file
fclose(fid);

% Rimuovi eventuali spazi extra
fileContent = strrep(fileContent, ' ', '');

% Suddividi il contenuto in righe
lines = strsplit(fileContent, ';');

% Inizializza una matrice vuota per i dati
data = [];

% Converti ogni riga in numeri
for i = 1:length(lines)
    if ~isempty(lines{i})  % Verifica che la riga non sia vuota
        numArray = str2double(strsplit(lines{i}, ','));
        data = [data; numArray];
    end
end

% Visualizza la matrice dei dati
disp(data);

%%
% NON VA
% Apri il file
fid = fopen('DecDec2.txt', 'rt');

% Inizializza una matrice vuota per i dati
data = [];

% Leggi ogni riga del file
while ~feof(fid)
    % Leggi la riga come stringa
    tline = fgetl(fid);
    if ischar(tline) && ~isempty(tline) % Assicurati che la riga non sia vuota
        % Rimuovi eventuali spazi extra
        tline = strrep(tline, ' ', '');
        % Dividi la riga in numeri separati da virgole
        numArray = str2double(strsplit(tline, ','));
        % Aggiungi la riga alla matrice dei dati
        data = [data; numArray];
    end
end

% Chiudi il file
fclose(fid);

% Visualizza la matrice dei dati
disp(data);

%%
% NON VA
% Apri il file
fid = fopen('DecDec2.txt', 'rt');

% Inizializza una matrice vuota per i dati
data = [];

% Leggi ogni riga del file
while ~feof(fid)
    % Leggi la riga come stringa
    tline = fgetl(fid);
    if ischar(tline) && ~isempty(tline) % Assicurati che la riga non sia vuota
        % Rimuovi eventuali spazi extra
        tline = strrep(tline, ' ', '');
        % Dividi la riga in numeri separati da virgole
        numArray = str2double(strsplit(tline, ','));
        % Verifica che la riga abbia esattamente 12 elementi
        if length(numArray) == 12
            % Aggiungi la riga alla matrice dei dati
            data = [data; numArray];
        else
            fprintf('Errore: Rilevata una riga con %d elementi anziché 12. Ignorata.\n', length(numArray));
        end
    end
end

% Chiudi il file
fclose(fid);

% Visualizza la matrice dei dati
disp(data);


