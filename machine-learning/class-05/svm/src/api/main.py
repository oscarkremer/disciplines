import cvxopt
import numpy as np
from sklearn.datasets.samples_generator import make_blobs
from sklearn.model_selection import train_test_split
from matplotlib import pyplot as plt
from sklearn.metrics import confusion_matrix
from src.models import SVM

if __name__=='__main__':
    svm = SVM()
    svm.fit(X_train, y_train)   
