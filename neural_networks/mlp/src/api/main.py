import random
import pandas as pd 
import numpy as np
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt

class Adaline:
    def __init__(self, learning_rate, sensibility=0.000001):
        self.learning_rate = learning_rate
        self.sensibility = sensibility
         
    def train(self, trainX, trainY):
        epochs = 1
        weights = np.random.rand(trainX.shape[1]+ 1)
        print('Pesos para Inicialização')
        print(weights)
        trainX = np.column_stack((-np.ones(trainX.shape[0]), trainX)) 
        predictions = np.zeros(trainY.shape[0])
        error_plot = []
        counter =1 
        while(True):
            old_error = self.error_mse(trainX, trainY, weights)
            for i in range(trainY.shape[0]):            
                predictions[i] = weights.dot(np.transpose(trainX))[i]
                weights+=self.learning_rate*(trainY[i] - predictions[i])*trainX[i]
            new_error = self.error_mse(trainX, trainY, weights)
            if  abs(old_error - new_error) <= self.sensibility:
                break
            error_plot.append(old_error)
            epochs+=1
        self.weights = weights
        self.epochs = epochs
        return error_plot

    def error_mse(self, trainX, trainY, weights):
        predictions = np.zeros(trainY.shape[0])
        for i in range(trainY.shape[0]):
            predictions[i] = weights.dot(np.transpose(trainX))[i]
        return np.sum(np.power(predictions-trainY, 2))/(predictions.shape[0])

    def predict(self, data):
        return self.activation(self.weights.dot(np.transpose(np.concatenate([[-1], data]))))

    def activation(self, x):
        return 1 if x >= 0  else -1


if __name__=='__main__':
    data = pd.read_csv('data.csv', decimal=',')
    train = data[:int(0.7*len(data))]
    test = data[int(0.7*len(data)):]
    trainX = np.array(train.drop(columns='d').values)
    trainY = np.array(train['d'].values)
    testX = np.array(test.drop(columns='d').values)    
    testY = np.array(test['d'].values)
    number_test = 5
    errors = []
    for i in range(number_test):
        model = Adaline(learning_rate=0.0025)
        errors.append(model.train(trainX, trainY))
        print('Pesos para após treinamento:')
        print(model.weights)
        print('Número total de épocas encontrado: {}'.format(model.epochs))
        predictions = []
        for data in testX:
            predictions.append(model.predict(np.array(data)))
        confusion = confusion_matrix(testY, predictions)
        print('Matriz de Confusao')
        print(confusion)
        validation = pd.read_csv('validation.csv').values
        for value in validation:
            print('--- Valores da Amostra')
            print(value)
            print('--- Predicao')
            print(model.predict(np.array(value)))
    
    for i in range(number_test):
       plt.plot(errors[i])
    plt.title('Mean Squared Error')
    plt.grid(True)
    plt.xlabel('Epochs')
    plt.ylabel('MSE')
    plt.savefig('plot.png')
    