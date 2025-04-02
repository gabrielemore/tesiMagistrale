function metriche_identificazione

% GoF media ML CIRC ON
GoF_tr = [26.4237 22.3797 34.0349 37.7199 37.3975 40.9629 48.5403 35.3713 27.5083 15.7198];
GoF_te = [24.521 23.3781 28.2513 37.5362 30.5862 45.4314 45.6975 34.5318 37.8828 18.5786];

mediana = median(GoF_tr);              % Calcola la mediana
percentili = prctile(GoF_tr, [25 75]); % Calcola il 25° e il 75° percentile

% Mostra i risultati
disp('----ML CIRC ON TRAIN----')
fprintf('Mediana: %.2f\n', mediana);
fprintf('25° percentile: %.2f\n', percentili(1));
fprintf('75° percentile: %.2f\n', percentili(2));


mediana = median(GoF_te);              % Calcola la mediana
percentili = prctile(GoF_te, [25 75]); % Calcola il 25° e il 75° percentile

% Mostra i risultati
disp('----ML CIRC ON TEST----')
fprintf('Mediana: %.2f\n', mediana);
fprintf('25° percentile: %.2f\n', percentili(1));
fprintf('75° percentile: %.2f\n', percentili(2));


% GoF media NL CIRC ON
GoF_tr = [52.163 49.0879 50.4733 55.1571 67.1608 56.7868 57.6734 60.6097 58.1176 52.3879];
GoF_te = [40.209 39.8906 37.186 57.3228 50.6695 58.5291 43.9576 55.0718 50.2657 54.9177];

mediana = median(GoF_tr);              % Calcola la mediana
percentili = prctile(GoF_tr, [25 75]); % Calcola il 25° e il 75° percentile

% Mostra i risultati
disp('----NL CIRC ON TRAIN----')
fprintf('Mediana: %.2f\n', mediana);
fprintf('25° percentile: %.2f\n', percentili(1));
fprintf('75° percentile: %.2f\n', percentili(2));


mediana = median(GoF_te);              % Calcola la mediana
percentili = prctile(GoF_te, [25 75]); % Calcola il 25° e il 75° percentile

% Mostra i risultati
disp('----NL CIRC ON TEST----')
fprintf('Mediana: %.2f\n', mediana);
fprintf('25° percentile: %.2f\n', percentili(1));
fprintf('75° percentile: %.2f\n', percentili(2));


% GoF media ML CIRC ON VAL. MESSORI
GoF_tr = [26.0636 0.068712 16.6339 30.0826 21.7021 28.0785 24.5171 26.2384 28.1374 10.7552];

mediana = median(GoF_tr);              % Calcola la mediana
percentili = prctile(GoF_tr, [25 75]); % Calcola il 25° e il 75° percentile

% Mostra i risultati
disp('----ML CIRC ON MESSORI----')
fprintf('Mediana: %.2f\n', mediana);
fprintf('25° percentile: %.2f\n', percentili(1));
fprintf('75° percentile: %.2f\n', percentili(2));


%GoF media NL CIRC ON VAL. MESSORI
GoF_tr = [40.7145 32.3272 37.4856 52.4835 33.0795 47.5924 35.0676 49.1995 43.9932 37.0929];

mediana = median(GoF_tr);              % Calcola la mediana
percentili = prctile(GoF_tr, [25 75]); % Calcola il 25° e il 75° percentile

% Mostra i risultati
disp('----NL CIRC ON MESSORI----')
fprintf('Mediana: %.2f\n', mediana);
fprintf('25° percentile: %.2f\n', percentili(1));
fprintf('75° percentile: %.2f\n', percentili(2));
end

