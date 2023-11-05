% EET-51 LAB 03
clear; clc;

% 1) Carregar as imagens deste experimento:
C = imread("baboon.tif");
B = imread("lena512.tif"); 
A = imread("peppers.tif");
% imshow(A);

% Calcule os histogramas para as três imagens
[histoA, binslocA] = imhist(A);
[histoB, binslocB] = imhist(B);
[histoC, binslocC] = imhist(C);

histA = histoA/(512^2);
histB = histoB/(512^2);
histC = histoC/(512^2);

% Plot dos histogramas com cores diferentes
hold on;
plot(binslocA, histA, 'r', 'LineWidth', 2);
plot(binslocB, histB, 'g', 'LineWidth', 2);
plot(binslocC, histC, 'b', 'LineWidth', 2);

% Personalize a área sob a curva com transparência
alpha = 0.3;
area(0:255, histA, 'FaceColor', 'r', 'FaceAlpha', alpha);
area(0:255, histB, 'FaceColor', 'g', 'FaceAlpha', alpha);
area(0:255, histC, 'FaceColor', 'b', 'FaceAlpha', alpha);

% Adicione legendas para as imagens
lgd = legend('baboon.tif', 'lena512.tif', 'peppers.tif', 'Location', 'Northwest');
set(lgd, 'FontSize', 12);

% Configurações de gráfico adicionais
title('Histograma', 'FontSize', 18);
xlabel('Valores dos pixels', 'FontSize', 18);
ylabel('Frequência relativa dos pixels', 'FontSize', 18);
xlim([0, 255]);

print('-depsc', 'hist.eps');

hold off; 

% 2) Cálculo da entropia de cada imagem:
imhist(A)
[counts_A,binLocations_A] = imhist(A);

Pr_A = counts_A/512^2;
Sum_A = sum(Pr_A)

plot(binLocations_A,Pr_A)
xlim([0 255]);

H_A = 0;

for i=1:256
    if histA(i) ~= 0
        H_A = H_A - histA(i)*log2(histA(i));
    end
end
H_A

H_B = 0;

for i=1:256
    if histB(i) ~= 0
        H_B = H_B - histB(i)*log2(histB(i));
    end
end
H_B

H_C = 0;

for i=1:256
    if histC(i) ~= 0
        H_C = H_C - histC(i)*log2(histC(i));
    end
end
H_C

% 3) Matrizes de erro de predição:
PA = zeros(512,512);
A = double(A);
for i=1:512
    for j=1:512
        if i == 1
            if j == 1
                PA(i,j) = A(i,j);
            else
                PA(i,j) = A(i,j) - A(i,j-1);
            end
        else
            if j == 1
                PA(i,j) = A(i,j) - A(i-1,j);
            else
                PA(i,j) = A(i,j) - A(i,j-1) - A(i-1,j) + A(i-1,j-1);
            end
        end
    end
end


h_PA = histogram(PA)
counts_PA = h_PA.Values
lims = h_PA.BinEdges;

Pr_PA = counts_PA/512^2;
Sum_PA = sum(Pr_PA)

plot(lims(1:end-1),Pr_PA)
xlim([lims(1) lims(end-1)]);

H_PA = 0;
for i=1:length(Pr_PA)
    if Pr_PA(i) ~= 0
        H_PA = H_PA - Pr_PA(i)*log2(Pr_PA(i));
    end
end
H_PA

% 4) Verificar que a imagem original pode ser obtida a partir da imagem transformada:
imshow(PA)

PA_rec = zeros(512,512);

for i=1:512
    for j=1:512
        if i == 1
            if j == 1
                PA_rec(i,j) = PA(i,j);
            else
                PA_rec(i,j) = PA(i,j) + PA_rec(i,j-1);
            end
        else
            if j == 1
                PA_rec(i,j) = PA(i,j) + PA_rec(i-1,j);
            else
                PA_rec(i,j) = PA(i,j) + PA_rec(i,j-1) + PA_rec(i-1,j) - PA_rec(i-1,j-1);
            end
        end
    end
end

imshow(PA_rec);
isequal(A, PA_rec)

% 5) Obtenção de uma nova imagem cujos pixels são os módulos dos pixels da imagem transformada:
PA2 = abs(PA);
imshow(PA2)

% 6) Nova entropia:
h_PA2 = histogram(PA2)
counts_PA2 = h_PA2.Values
lims2 = h_PA2.BinEdges;

Pr_PA2 = counts_PA2/512^2;
Sum_PA2 = sum(Pr_PA2)

plot(lims2(1:end-1),Pr_PA2)
xlim([lims2(1) lims2(end-1)]);

H_PA2 = 0;

for i=1:length(Pr_PA2)
    if Pr_PA2(i) ~= 0
        H_PA2 = H_PA2 - Pr_PA2(i)*log2(Pr_PA2(i));
    end
end
H_PA2 % Sem contar o bit de sinal!

PA3 = sign(PA);
h_PA3 = histogram(PA3)
counts_PA3 = h_PA3.Values
lims3 = h_PA3.BinEdges;

Pr_PA3 = counts_PA3/sum(counts_PA3);
Sum_PA3 = sum(Pr_PA3)

plot(lims3(1:end-1),Pr_PA3)
xlim([lims3(1) lims3(end-1)]);

H_PA3 = 0;

s1 = Pr_PA3(1)
s2 = Pr_PA3(2) + Pr_PA3(3)
H_PA3_s = - s1*log2(s1)  - s2*log2(s2) 
for i=1:length(Pr_PA3)
    if Pr_PA3(i) ~= 0
        H_PA3 = H_PA3 - Pr_PA3(i)*log2(Pr_PA3(i));
    end
end
H_PA3 % Sem contar o bit de sinal!

% 7) Codificação de Golomb:
M = sum(PA2,'all')/512^2
m = ceil(log2(M/2))

GA = zeros(512,512);
for i=1:512
    for j=1:512
         GA(i,j) = 2 + m + floor(PA2(i,j)/2^m);
    end
end

bits_sum = sum(GA,'all')
reduction = 100*(1 - bits_sum/(8*512^2))

function golomb_encoded = golomb_encode(values, m)
    golomb_encoded = "";
    
    for i = 1:length(values)
        value = values(i);
        
        if value < 0
            error('O valor deve ser não negativo.');
        end
        
        q = floor(value / (2^m));
        r = value - q * (2^m);
        
        unary_part = repmat('0', 1, q);
        binary_part = dec2bin(r, m);
        
        golomb_encoded = strcat(golomb_encoded, unary_part, '1', binary_part);
    end
end
