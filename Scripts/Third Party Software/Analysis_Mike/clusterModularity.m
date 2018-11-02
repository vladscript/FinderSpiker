function [Ci]=clusterModularity(matriz_adyacente,nsim)
    matriz_conexiones=zeros(size(matriz_adyacente));
    num_neuronas = length(matriz_adyacente);
    for i=1:nsim
       [Ci,Q]=modularity_und(matriz_adyacente);
        for j=1:num_neuronas
            modulo = Ci(j);
            indices = find(Ci==modulo);
            matriz_conexiones(j,indices)=matriz_conexiones(j,indices)+1;
        end
        disp(i/nsim)
    end
    for i=1:num_neuronas
       matriz_conexiones(i,i)=0;
    end
    [Ci,Q]=modularity_und(matriz_conexiones);
%     figure(2);imagesc(matriz_conexiones);colorbar;
end
