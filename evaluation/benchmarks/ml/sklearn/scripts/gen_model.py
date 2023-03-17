from sklearn.linear_model import LogisticRegression
import pickle
import sys

reg = LogisticRegression(max_iter=int(sys.argv[1]), 
                         solver='newton-cholesky', 
                         multi_class='ovr')
with open('./tmp/model.obj', 'w+b') as file:
    pickle.dump(reg, file)