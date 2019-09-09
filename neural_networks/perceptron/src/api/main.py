import random
import pandas as pd 
import numpy as np
from sklearn.metrics import confusion_matrix

class Perceptron:
    def __init__(self, learning_rate):
        self.learning_rate = learning_rate
        
    def train(self, trainX, trainY):
        epochs = 1
        weights = np.random.rand(trainX.shape[1]+ 1)
        print('Pesos para Inicialização')
        print(weights)
        trainX = np.column_stack((-np.ones(trainX.shape[0]), trainX)) 
        predictions = np.zeros(trainY.shape[0])
        while(True): 
            for i in range(trainY.shape[0]):
                predictions[i] = self.activation(weights.dot(np.transpose(trainX))[i])
                if predictions[i] != trainY[i]:
                    weights+=self.learning_rate*(trainY[i] - predictions[i])*trainX[i]
            if np.array_equal(predictions, trainY):
                break
            epochs+=1
        self.weights = weights
        self.epochs = epochs

    def predict(self, data):
        return self.activation(self.weights.dot(np.transpose(np.concatenate([[-1], data]))))

    def activation(self, x):
        return 1 if x >= 0  else -1


if __name__=='__main__':
    data = pd.read_csv('data.csv', delimiter=';')
    data = data.iloc[np.random.permutation(len(data))]
    train = data[:int(0.7*len(data))]
    test = data[int(0.7*len(data)):]
    trainX = np.array(train.drop(columns='d').values)
    trainY = np.array(train['d'].values)
    testX = np.array(test.drop(columns='d').values)
    testY = np.array(test['d'].values)
    model = Perceptron(0.01)
    model.train(trainX, trainY)
    print('Pesos para após treinamento:')
    print(model.weights)
    print('Número total de épocas encontrado:')
    print(model.epochs)
    predictions = []
    for data in testX:
        predictions.append(model.predict(np.array(data)))
    confusion = confusion_matrix(testY, predictions)
    print('Matriz de Confusao')
    print(confusion)
    validation = pd.read_csv('validation.csv').values
    for value in validation:
        print('valores da amostra')
        print(value)
        print('predicao')
        print(model.predict(np.array(value)))
