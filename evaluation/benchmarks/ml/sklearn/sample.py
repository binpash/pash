from sklearn.model_selection import train_test_split
from sklearn import datasets
from sklearn.linear_model import LogisticRegression

max_iter=10

dataset = datasets.fetch_covtype()
X_train, X_test, y_train, y_test = train_test_split(dataset.data, dataset.target, 
                                                    test_size=0.2, 
                                                    random_state=0)
control_model = LogisticRegression(max_iter=max_iter, 
                                   solver='newton-cholesky', 
                                   multi_class='ovr')
control_model.fit(X_train, y_train)
control_score = control_model.score(X_test, y_test)
print(control_score)