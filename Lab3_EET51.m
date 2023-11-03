%% EET-51 LAB 03
clear; clc;

%% 1) Carregar as imagens deste experimento:
A = imread("baboon.tif");
% A = imread("lena512.tif"); 
% A = imread("peppers.tif");
imshow(A);

%% 2) Cálculo da entropia de cada imagem:
imhist(A)
[counts_A,binLocations_A] = imhist(A);

Pr_A = counts_A/512^2;
Sum_A = sum(Pr_A)

plot(binLocations_A,Pr_A)
xlim([0 255]);

H_A = 0;

for i=1:256
    if Pr_A(i) ~= 0
        H_A = H_A - Pr_A(i)*log2(Pr_A(i));
    end
end
H_A

%% 3) Matrizes de erro de predição:
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

for i=1:314
    if Pr_PA(i) ~= 0
        H_PA = H_PA - Pr_PA(i)*log2(Pr_PA(i));
    end
end
H_PA

%% 4) Verificar que a imagem original pode ser obtida a partir da imagem transformada:
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

isequal(A, PA_rec)

%% 5) Obtenção de uma nova imagem cujos pixels são os módulos dos pixels da imagem transformada:
PA2 = abs(PA);
imshow(PA2)

%% 6) Nova entropia:
h_PA2 = histogram(PA2)
counts_PA2 = h_PA2.Values
lims2 = h_PA2.BinEdges;

Pr_PA2 = counts_PA2/512^2;
Sum_PA2 = sum(Pr_PA2)

plot(lims2(1:end-1),Pr_PA2)
xlim([lims2(1) lims2(end-1)]);

H_PA2 = 0;

for i=1:166
    if Pr_PA2(i) ~= 0
        H_PA2 = H_PA2 - Pr_PA2(i)*log2(Pr_PA2(i));
    end
end
H_PA2 % Sem contar o bit de sinal!

%% 7) Codificação de Golomb:
M = sum(PA2,'all')/512^2;
m = ceil(log2(M/2))

GA = zeros(512,512);
for i=1:512
    for j=1:512
         GA(i,j) = 2 + m + floor(PA2(i,j)/2^m);
    end
end

bits_sum = sum(GA,'all')
reduction = 100*(1 - bits_sum/(8*512^2))
