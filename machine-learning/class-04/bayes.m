%%%% Código para disciplina de aprendizado de máquina
%%%% Oscar Schmitt Kremer - Professor : Ms Lucian Schiavon
%%%% Metodo de Naive Bayes
%% Limpar console, variaveis e fechar janelas
clc
clear all;
close all;
%% Carregar e separar dados de treino e teste
data = load('data.mat');
trainX = data.X;
trainY = data.Y;
testX = data.newX;
testY = data.newY;
%% Inicializar alguns vetores como zero e constantes utilizadas
count = zeros(5,2);
countword = zeros(5,1);
message = zeros(2,1);
num_word = 20;
train_length = size(trainX);train_length = train_length(1);
test_length = size(testX); test_length = test_length(1);
total_words = num_word*train_length;
%% Contagem de spans e ocorrencia das palavaras em spans ou hams
for i=1:train_length
    if trainY(i)==-1
        message(1) = message(1)+1;
    else
        message(2) = message(2)+1;
    end
    for j=1:5
        if trainX(i,j)~=0
            countword(j) = countword(j) + trainX(i,j);
            if trainY(i)==1
                count(j,2) = count(j,2) + trainX(i,j);
            else
                count(j,1) = count(j,1) + trainX(i,j);        
            end
        end
    end
end
%% Calculo das probabilidades a partir dos valores contados
ham_or_spam = message/sum(message);
count = count/total_words;
w_ham_spam = [count(:,1)/ham_or_spam(1) count(:,2)/ham_or_spam(2)];
word = countword/total_words;
result = zeros(test_length);
%% Laço final para classificar emails de test a partir da análise feita
for i=1:test_length
    p_message = prob(word, testX(i,:)); 
    p_ham_e = ham_or_spam(1)*prob(w_ham_spam(:,1), testX(i,:))/p_message; 
    p_spam_e = ham_or_spam(2)*prob(w_ham_spam(:,2), testX(i,:))/p_message;
    if p_spam_e > p_ham_e
        result(i) =  1;
    else
        result(i) = -1;
    end 
end
%% Encontrar acertos na forma de uma matriz de confusão
final = evaluate(result, testY);
accuracy = (final(1,1)+final(2,2))/test_length
true_positive_rate = (final(2,2))/(final(2,2)+final(2,1))
true_false_rate = (final(1,1))/(final(1,1)+final(1,2))
