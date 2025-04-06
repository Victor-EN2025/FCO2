function ax = ptable(rows, cols)
    % Crée une figure avec un arrangement personnalisé de sous-graphes
    % rows et cols sont des vecteurs spécifiant la disposition des sous-figures
    ax = zeros(rows(1), rows(2)); % Initialiser un tableau pour les axes
    for i = 1:rows(1)
        for j = 1:cols(i)
            ax(i,j) = subplot(rows(1), cols(i), (i-1)*cols(i) + j);
        end
    end
end
