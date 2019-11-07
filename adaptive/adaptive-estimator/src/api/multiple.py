import os
import numpy as np 
import matplotlib.pyplot as plt 
from src.models import Estimator
from random import *

def f(t):
    sig = []
    for i in range(t.shape[0]):
        if i<t.shape[0]/3:
            sig.append(1)
        else:
            if i< 2*t.shape[0]/3:
                sig.append(0)
            else:
                sig.append(1)
    return np.array(sig[:t.shape[0]])


def calculate_terms(array_model, delta_time):
    print(array_model.shape)   
    a = np.zeros(array_model.shape[1])
    b = np.zeros(array_model.shape[1])
    for i in range(array_model.shape[1]):
        b[i] = (array_model[0][i]-1)/delta_time
        a[i] = (array_model[1][i])/array_model[0][i]
    return a,b
#    b = (array_model[0]-1)/delta_time
#    a = (array_model[1])
#    return np.array(sig[:t.shape[0]])


if __name__=='__main__':
    u, y ,models, thetas, p = [], [], [], [], []
    step = 0.001
    final_time = 5.0
    t = np.arange(0, step + final_time, step)
    update_times = [5.1]
    a = 1
    b = 2
    sig = f(t)

    for i in range(t.shape[0]):
        if i == 0:
            y.append(0)
            u.append([0, 0])
        else:
            y.append((1/(1+b*step))*y[i-1] + (a/(1+b*step))*sig[i])        
            u.append([y[i-1], sig[i]])
    u = np.transpose(np.array(u))
    analysis = 'error'
#    analysis = 'mean'
    y = np.array(y)
    for update_time in update_times:
        model = Estimator(update_time, step + final_time)
        model.train(t, u, y)
        models.append(model)

    for i in range(t.shape[0]):
        error_list  = []
        if analysis=='error':
            for model in models:
                error_list.append(model.error[i])
            index = np.argmin(np.array(error_list))
        else:
            for model in models:
                error_list.append(model.mean[i])
            index = np.argmin(np.array(error_list))
   #     print(models[index])

        thetas.append(models[index].theta_plot[i])
        
        p.append(models[index].p[i])

    y_hat  = np.zeros(t.shape[0])
    error = np.zeros(t.shape[0])
    for i in range(t.shape[0]):
        y_hat[i] = thetas[i][0]*u[0][i] + thetas[i][1]*u[1][i]
        error[i] = (y[i] - y_hat[i])/y[i]

    plt.plot(t, y)
    for model in models:
        plt.plot(t, model.y_hat)
    plt.show()

    for model in models:
        plt.plot(t, np.transpose(np.array(model.theta_plot))[0])
        plt.plot(t, np.transpose(np.array(model.theta_plot))[1])
        plt.show()

    a, b = calculate_terms(np.transpose(np.array(model.theta_plot)), step)
    plt.plot(a)
    plt.show()

    plt.plot(b)
    plt.show()
    