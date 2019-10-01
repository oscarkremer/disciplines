import scipy.io
import numpy as np
from sklearn.preprocessing import StandardScaler
from matplotlib import pyplot as plt

class SVM:
  
    def __init__(self, X, y, C, kernel, alphas, b, errors):
        self.X = X               
        self.y = y               
        self.C = C               
        self.kernel = kernel     
        self.alphas = alphas     
        self.b = b               
        self.errors = errors     
        self._obj = []           
        self.m = len(self.X)     

def menu():
    print('[1] - Kernel Polinomial de Grau 1 \n')
    print('[2] - Kernel Polinomial de Grau 2 \n')
    print('[3] - Kernel Polinomial de Grau 3 \n')
    print('[Outros] - Kernel Gaussiano \n')
    kernel = input('Entre com sua escolha: ')
    try:
        if int(kernel) == 1:
            return 'Linear'
        else:
            if int(kernel) == 2:
                return 'Quadratic'
            else:
                if int(kernel) == 3:
                    return 'Cubic'
                else:
                    return 'Gaussian'
    except:
        return 'Gaussian'

def linear_kernel(x, y):
    return x @ y.T

def quadratic_kernel(x, y):
    kernel = x @ y.T
    if np.ndim(x) == 1 and np.ndim(y) == 1:
        result = kernel**2
    elif (np.ndim(x) > 1 and np.ndim(y) == 1) or (np.ndim(x) == 1 and np.ndim(y) > 1):
        result = np.power(kernel,2)
    elif np.ndim(x) > 1 and np.ndim(y) > 1:
        result = kernel.dot(kernel)
    return result

def cubic_kernel(x, y):
    kernel = x @ y.T
    if np.ndim(x) == 1 and np.ndim(y) == 1:
        result = kernel**3
    elif (np.ndim(x) > 1 and np.ndim(y) == 1) or (np.ndim(x) == 1 and np.ndim(y) > 1):
        result = np.power(kernel,3)
    elif np.ndim(x) > 1 and np.ndim(y) > 1:
        result = (kernel.dot(kernel)).dot(kernel)
    return result

def gaussian_kernel(x, y, sigma=1):
    if np.ndim(x) == 1 and np.ndim(y) == 1:
        result = np.exp(- (np.linalg.norm(x - y, 2)) ** 2 / (2 * sigma ** 2))
    elif (np.ndim(x) > 1 and np.ndim(y) == 1) or (np.ndim(x) == 1 and np.ndim(y) > 1):
        result = np.exp(- (np.linalg.norm(x - y, 2, axis=1) ** 2) / (2 * sigma ** 2))
    elif np.ndim(x) > 1 and np.ndim(y) > 1:
        result = np.exp(- (np.linalg.norm(x[:, np.newaxis] - y[np.newaxis, :], 2, axis=2) ** 2) / (2 * sigma ** 2))
    return result

def objective_function(alphas, target, kernel, X_train):
    return np.sum(alphas)-0.5*np.sum((target[:, None]*target[None, :])*kernel(X_train, X_train)*(alphas[:, None]*alphas[None, :]))

def decision_function(alphas, target, kernel, X_train, x_test, b):
    return (alphas * target) @ kernel(X_train, x_test) - b
    
def plot_decision_boundary(model, ax, type_kernel, resolution=400, colors=('b', 'k', 'r'), levels=(-1, 0, 1)):
        ax[0].scatter(X_train[:, 0], X_train[:, 1], c=y, cmap=plt.cm.viridis)
        xrange = np.linspace(model.X[:,0].min(), model.X[:,0].max(), resolution)
        yrange = np.linspace(model.X[:,1].min(), model.X[:,1].max(), resolution)
        grid = [[decision_function(model.alphas, model.y,
                                   model.kernel, model.X,
                                   np.array([xr, yr]), model.b) for xr in xrange] for yr in yrange]
        grid = np.array(grid).reshape(len(xrange), len(yrange))
        ax[0].grid(True)
        ax[0].set_title('Training Data')
        ax[1].grid(True)
        ax[1].contour(xrange, yrange, grid, levels=levels, linewidths=(1, 1, 1),
                   linestyles=('--', '-', '--'), colors=colors)
        ax[1].scatter(model.X[:,0], model.X[:,1],
                   c=model.y, cmap=plt.cm.viridis, lw=0, alpha=0.75)        
        mask = np.round(model.alphas, decimals=2) != 0.0
        ax[1].scatter(model.X[mask,0], model.X[mask,1],
                   c=model.y[mask], cmap=plt.cm.viridis, lw=1, edgecolors='k')
        ax[1].set_title('Classification - SVM - Kernel - {}'.format(type_kernel))
        return grid, ax

def take_step(i1, i2, model):
    if i1 == i2:
        return 0, model
    alph1 = model.alphas[i1]
    alph2 = model.alphas[i2]
    y1 = model.y[i1]
    y2 = model.y[i2]
    E1 = model.errors[i1]
    E2 = model.errors[i2]
    s = y1 * y2
    if (y1 != y2):
        L = max(0, alph2 - alph1)
        H = min(model.C, model.C + alph2 - alph1)
    elif (y1 == y2):
        L = max(0, alph1 + alph2 - model.C)
        H = min(model.C, alph1 + alph2)
    if (L == H):
        return 0, model
    k11 = model.kernel(model.X[i1], model.X[i1])
    k12 = model.kernel(model.X[i1], model.X[i2])
    k22 = model.kernel(model.X[i2], model.X[i2])
    eta = 2*k12 - k11 - k22
    if (eta < 0):
        a2 = alph2 - y2 * (E1 - E2) / eta
        if L < a2 < H:
            a2 = a2
        elif (a2 <= L):
            a2 = L
        elif (a2 >= H):
            a2 = H
    else:
        alphas_adj = model.alphas.copy()
        alphas_adj[i2] = L
        Lobj = objective_function(alphas_adj, model.y, model.kernel, model.X) 
        alphas_adj[i2] = H
        Hobj = objective_function(alphas_adj, model.y, model.kernel, model.X)
        if Lobj > (Hobj + eps):
            a2 = L
        elif Lobj < (Hobj - eps):
            a2 = H
        else:
            a2 = alph2
    if a2 < 1e-8:
        a2 = 0.0
    elif a2 > (model.C - 1e-8):
        a2 = model.C    
    if (np.abs(a2 - alph2) < eps * (a2 + alph2 + eps)):
        return 0, model
    a1 = alph1 + s * (alph2 - a2)
    b1 = E1 + y1 * (a1 - alph1) * k11 + y2 * (a2 - alph2) * k12 + model.b
    b2 = E2 + y1 * (a1 - alph1) * k12 + y2 * (a2 - alph2) * k22 + model.b
    if 0 < a1 and a1 < C:
        b_new = b1
    elif 0 < a2 and a2 < C:
        b_new = b2
    else:
        b_new = (b1 + b2)*0.5
    model.alphas[i1] = a1
    model.alphas[i2] = a2
    for index, alph in zip([i1, i2], [a1, a2]):
        if 0.0 < alph < model.C:
            model.errors[index] = 0.0    
    non_opt = [n for n in range(model.m) if (n != i1 and n != i2)]
    model.errors[non_opt] = model.errors[non_opt] + \
                            y1*(a1 - alph1)*model.kernel(model.X[i1], model.X[non_opt]) + \
                            y2*(a2 - alph2)*model.kernel(model.X[i2], model.X[non_opt]) + model.b - b_new
    
    model.b = b_new
    return 1, model

def examine_example(i2, model):
    y2 = model.y[i2]
    alph2 = model.alphas[i2]
    E2 = model.errors[i2]
    r2 = E2 * y2
    if ((r2 < -tol and alph2 < model.C) or (r2 > tol and alph2 > 0)):
        if len(model.alphas[(model.alphas != 0) & (model.alphas != model.C)]) > 1:
            if model.errors[i2] > 0:
                i1 = np.argmin(model.errors)
            elif model.errors[i2] <= 0:
                i1 = np.argmax(model.errors)
            step_result, model = take_step(i1, i2, model)
            if step_result:
                return 1, model
        for i1 in np.roll(np.where((model.alphas != 0) & (model.alphas != model.C))[0],
                          np.random.choice(np.arange(model.m))):
            step_result, model = take_step(i1, i2, model)
            if step_result:
                return 1, model
        for i1 in np.roll(np.arange(model.m), np.random.choice(np.arange(model.m))):
            step_result, model = take_step(i1, i2, model)
            if step_result:
                return 1, model
    return 0, model

def train(model):
    numChanged = 0
    examineAll = 1
    while(numChanged > 0) or (examineAll):
        numChanged = 0
        if examineAll:
            for i in range(model.alphas.shape[0]):
                examine_result, model = examine_example(i, model)
                numChanged += examine_result
                if examine_result:
                    obj_result = objective_function(model.alphas, model.y, model.kernel, model.X)
                    model._obj.append(obj_result)
        else:
            for i in np.where((model.alphas != 0) & (model.alphas != model.C))[0]:
                examine_result, model = examine_example(i, model)
                numChanged += examine_result
                if examine_result:
                    obj_result = objective_function(model.alphas, model.y, model.kernel, model.X)
                    model._obj.append(obj_result)
        if examineAll == 1:
            examineAll = 0
        elif numChanged == 0:
            examineAll = 1
    return model

if __name__=='__main__':
    mat = scipy.io.loadmat('data/raw/data.mat')
    data = mat['data']
    X_train = data[:,:2]
    y = data[:,2]
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train, y)
    C = 0.5
    initial_b = 0.0
    tol = 0.01
    eps = 0.01
    m = len(X_train_scaled)
    initial_alphas = np.zeros(m)
    type_kernel = menu()
    if type_kernel=='Linear':
        model = SVM(X_train_scaled, y, C, linear_kernel,
                    initial_alphas, initial_b, np.zeros(m))
    if type_kernel=='Quadratic':
        model = SVM(X_train_scaled, y, C, quadratic_kernel,
                    initial_alphas, initial_b, np.zeros(m))
    if type_kernel=='Cubic':
        model = SVM(X_train_scaled, y, C, cubic_kernel,
                    initial_alphas, initial_b, np.zeros(m))
    if type_kernel=='Gaussian':
        model = SVM(X_train_scaled, y, C, gaussian_kernel,
                    initial_alphas, initial_b, np.zeros(m))
    
    initial_error = decision_function(model.alphas, model.y, model.kernel,
                                    model.X, model.X, model.b) - model.y
    model.errors = initial_error
    np.random.seed(0)
    output = train(model)
    fig, ax = plt.subplots(2)
    grid, ax = plot_decision_boundary(output, ax, type_kernel)
    plt.show()
