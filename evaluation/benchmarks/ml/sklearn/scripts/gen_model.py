from sklearn.linear_model import LogisticRegression
import pickle

reg = LogisticRegression(max_iter=10000)
with open('./tmp/model.obj', 'w+b') as file:
    pickle.dump(reg, file)