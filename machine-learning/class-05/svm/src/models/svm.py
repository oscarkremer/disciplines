
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
