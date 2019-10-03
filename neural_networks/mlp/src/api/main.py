import random
import pandas as pd 
import numpy as np
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt

class MLP:
    def __init__(self, number_neurons, learning_rate, momentum, sensibility=0.000001):
        self.number_neurons = number_neurons
        self.sensibility = sensibility
        self.alpha = momentum
        self.eta = learning_rate
         
    def train(self, trainX, trainY):
        epochs = 1
        weights_1 = np.random.rand(self.number_neurons[1], self.number_neurons[0], 3)
        weights_2 = np.random.rand(self.number_neurons[2], self.number_neurons[1], 3)
        l_1 = np.zeros(weights_1.shape[0])
        l_2 = np.zeros(weights_2.shape[0])
        Y_1 = np.zeros(weights_1.shape[0])
        Y_2 = np.zeros(weights_2.shape[0])
        for k in range(trainY.shape[0]):            
            for j in range(l_1.shape[0]):
                l_1[j] = self.sum_weights(weights_1[j,:,1], trainX[k])
                Y_1[j] = self.activation(l_1[j])
            for j in range(l_2.shape[0]):
                l_2[j] = self.sum_weights(weights_2[j,:,1], trainX[k])
                Y_2[j] = self.activation(l_2[j])
        old_error = self.error_mse(trainX, trainY, weights_1, weights_2)
        new_error = old_error
        while(True):
            old_error = new_error
            for k in range(trainY.shape[0]):            
                for j in range(l_1.shape[0]):
                    l_1[j] = self.sum_weights(weights_1[j,:,1], trainX[k])
                    Y_1[j] = self.activation(l_1[j]) 
                for j in range(l_2.shape[0]):
                    l_2[j] = self.sum_weights(weights_2[j,:,1], Y_1)
                    Y_2[j] = self.activation(l_2[j]) 

                delta_1 = np.zeros(weights_1.shape[0])
                delta_2 = np.zeros(weights_2.shape[0])
                for j in range(weights_2.shape[0]):
                    delta_2[j] = (trainY[k] - Y_2[j])*self.dev_activation(l_2[j])
                    for i in range(weights_2.shape[1]):
                        weights_2[j][i][2] = (self.alpha+1)*weights_2[j][i][1] - self.alpha*weights_2[j][i][0] + self.eta*delta_2[j]*Y_1[i]
                        weights_2[j][i][0] = weights_2[j][i][1] 
                        weights_2[j][i][1] = weights_2[j][i][2]     
         
                for j in range(weights_1.shape[0]):
                    delta_1[j] = self.calculate_sum(delta_2, weights_2[:,j,1])*self.dev_activation(l_1[j])
                    for i in range(weights_1.shape[1]):
                        weights_1[j][i][2] = (self.alpha+1)*weights_1[j][i][1] - self.alpha*weights_1[j][i][0] + self.eta*delta_1[j]*trainX[k][i]
                        weights_1[j][i][0] = weights_1[j][i][1]             
                        weights_1[j][i][1] = weights_1[j][i][2]                         
                    
            new_error = self.error_mse(trainX, trainY, weights_1, weights_2)
            
            print('Error - {}'.format(new_error))
            print(epochs)
            epochs+=1
            if abs(new_error - old_error) < self.sensibility:
                break
        self.weights_1 = weights_1
        self.weights_2 = weights_2
        

    def error_mse(self, trainX, trainY, weights_1, weights_2):
        prediction = np.zeros(trainX.shape[0])
        l_1 = np.zeros(weights_1.shape[0])
        l_2 = np.zeros(weights_2.shape[0])
        Y_1 = np.zeros(weights_1.shape[0])
        Y_2 = np.zeros(weights_2.shape[0])
        for i in range(trainX.shape[0]):
            for j in range(l_1.shape[0]):
                l_1[j] = self.sum_weights(weights_1[j,:,1], trainX[i]) 
                Y_1[j] = self.activation(l_1[j])
            for j in range(l_2.shape[0]):
                l_2[j] = self.sum_weights(weights_2[j,:,1], Y_1)
                Y_2[j] = self.activation(l_2[j])
            prediction[i] = Y_2
        return 0.5*np.sum(np.power(prediction-trainY, 2))/(prediction.shape[0])

    def predict(self, data):
        l_1 = np.zeros(self.weights_1.shape[0])
        l_2 = np.zeros(self.weights_2.shape[0])
        Y_1 = np.zeros(self.weights_1.shape[0])
        Y_2 = np.zeros(self.weights_2.shape[0])
        for j in range(l_1.shape[0]):
            l_1[j] = self.sum_weights(self.weights_1[j,:,1], data) 
            Y_1[j] = self.activation(l_1[j])
        for j in range(l_2.shape[0]):
            l_2[j] = self.sum_weights(self.weights_2[j,:,1], data) 
            Y_2[j] = self.activation(l_2[j])
        return Y_2
        
    def activation(self, x):
        return 1/(1+np.exp(-x))

    def dev_activation(self, x):
        return np.exp(-x)/((1+np.exp(-x))**2)

    def calculate_sum(self, deltas, weights):
        sum_ = 0
        for i in range(deltas.shape[0]):
            sum_+=deltas[i]*weights[i]
        return sum_

    def sum_weights(self, weights, trainX):
        sum_ = 0
        for i in range(trainX.shape[0]):
            sum_ += trainX[i]*weights[i]
        return sum_

if __name__=='__main__':
    train = pd.read_csv('data/train/project_1.csv', decimal=',')
    trainX = np.array(train.drop(columns='d').values)
    trainY = np.array(train['d'].values)
    model = MLP([3, 10, 1], 0.1, 0.8, sensibility=0.000001)
    model.train(trainX, trainY)
    for i in range(trainY.shape[0]):
        print('{0} - {1}'.format(trainY[i], model.predict(trainX[i])))
#    print('Número total de épocas encontrado: {}'.format(model.epochs))
#    predictions = []
#    for data in testX:
#        predictions.append(model.predict(np.array(data)))
#    confusion = confusion_matrix(testY, predictions)
#    print('Matriz de Confusao')
#    print(confusion)
#    validation = pd.read_csv('validation.csv').values
#    for value in validation:
#        print('--- Valores da Amostra')
#        print(value)
#        print('--- Predicao')
#        print(model.predict(np.array(value)))
    
    