import random
import pandas as pd 
import numpy as np
from sklearn.metrics import mean_squared_error 
import matplotlib.pyplot as plt

        
class MLP:
    def __init__(self, number_neurons, learning_rate, momentum, sensibility=0.000001):
        self.number_neurons = number_neurons
        self.sensibility = sensibility
        self.alpha = momentum
        self.eta = learning_rate
        self.initialize_weights()

    def initialize_weights(self):
        self.weights, self.bias = [],[]
        self.weights.append(np.random.rand(self.number_neurons[1], self.number_neurons[0], 3))
        self.weights.append(np.random.rand(self.number_neurons[2], self.number_neurons[1], 3))
        self.bias.append(np.random.rand(self.number_neurons[1], 3))
        self.bias.append(np.random.rand(self.number_neurons[2], 3))
        self.weights = np.array(self.weights)            
        self.bias = np.array(self.bias)        

    def fit(self, X, Y):
        self.epochs = 1
        pred = np.zeros(Y.shape[0])
        for j in range(Y.shape[0]):            
            pred[j] = self.predict(X[j])
        self.error = []
        old_error = mean_squared_error(Y, pred)
        self.error.append(old_error)
        new_error = old_error
        while(True):
            old_error = new_error
            for j in range(Y.shape[0]):            
                l_1 = self.weights[0][:,:,1].dot(np.transpose(X[j])) - self.bias[0][:,1]
                Y_1 = self.activation(l_1) 
                l_2 = self.weights[1][:,:,1].dot(np.transpose(Y_1))  - self.bias[1][:,1]
                Y_2 = self.activation(l_2) 
                delta = (Y[j] - Y_2)*self.dev_activation(l_2)
                self.bias[1][:,2] = (self.alpha+1)*self.bias[1][:,1] - self.alpha*self.bias[1][:,0] - self.eta*delta
                self.weights[1][:,:,2] = (self.alpha+1)*self.weights[1][:,:,1] - self.alpha*self.weights[1][:,:,0] + self.eta*delta*Y_1 
                delta = np.transpose(self.weights[1][:,:,1]).dot((delta))*self.dev_activation(l_1)
                self.bias[0][:,2] = (self.alpha+1)*self.bias[0][:,1] - self.alpha*self.bias[0][:,0] - self.eta*delta
                self.weights[0][:,:,2] = (self.alpha+1)*self.weights[0][:,:,1] - self.alpha*self.weights[0][:,:,0] + self.eta*(delta.reshape(delta.shape[0],1)).dot(np.transpose(X[j,:].reshape(X[j,:].shape[0],1)))                
                for i in range(2):
                    self.bias[i][:,0] = self.bias[i][:,1] 
                    self.bias[i][:,1] = self.bias[i][:,2]     
                    self.weights[i][:,:,0] = self.weights[i][:,:,1]             
                    self.weights[i][:,:,1] = self.weights[i][:,:,2]                         
             
            for j in range(Y.shape[0]):            
                pred[j] = self.predict(X[j])      
            new_error = mean_squared_error(Y, pred)
            self.error.append(new_error)
            self.epochs+=1
            print(new_error)
        
            if abs(new_error - old_error) < self.sensibility:
                print(new_error)
                break

    def predict(self, data):
        Y_1 = self.activation(self.weights[0][:,:,1].dot(np.transpose(data)) - self.bias[0][:,1])
        return self.activation(self.weights[1][:,:,1].dot(np.transpose(Y_1)) - self.bias[1][:,1])
        
    def activation(self, x):
        return 1/(1+np.exp(-x))

    def dev_activation(self, x):
        return np.exp(-x)/((1+np.exp(-x))**2)

if __name__=='__main__':
    train = pd.read_csv('data/train/project_1.csv', decimal=',')
    trainX = np.array(train.drop(columns='d').values)
    trainY = np.array(train['d'].values)    
    test = pd.read_csv('data/test/project_1.csv', decimal=',')
    testX = np.array(test.drop(columns='d').values)
    testY = np.array(test['d'].values)    

    model = MLP([3, 10, 1], 0.1, 0.8, sensibility=0.000001)
    model.fit(trainX, trainY)
    pred = []
    print('------ TESTE--------')
    for i in range(testY.shape[0]):
        pred.append(model.predict(testX[i]))
    pred = np.array(pred).reshape((testY.shape[0]))
    print('EPOCAS - {}'.format(model.epochs))    
    error_test = (testY - pred)
    print('ERROR MEDIO TESTE - {}'.format(np.mean(error_test)))
    print('VARIANCIA TESTE - {}'.format(np.std(error_test)**2))
    plt.plot(model.error)
    plt.title('Mean Squared Error')
    plt.grid(True)
    plt.xlabel('Epochs')
    plt.ylabel('MSE')
    plt.savefig('plot.png')
    