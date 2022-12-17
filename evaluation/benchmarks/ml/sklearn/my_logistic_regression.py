from sklearn.linear_model import _logistic

import warnings
import numpy as np

class MyLogisticRegression(_logistic.LogisticRegression):
    def warning_checks(self):
        if not isinstance(self.C, _logistic.numbers.Number) or self.C < 0:
            raise ValueError("Penalty term must be positive; got (C=%r)" % self.C)
        if self.penalty == "elasticnet":
            if (
                not isinstance(self.l1_ratio, _logistic.numbers.Number)
                or self.l1_ratio < 0
                or self.l1_ratio > 1
            ):
                raise ValueError(
                    "l1_ratio must be between 0 and 1; got (l1_ratio=%r)"
                    % self.l1_ratio
                )
        elif self.l1_ratio is not None:
            warnings.warn(
                "l1_ratio parameter is only used when penalty is "
                "'elasticnet'. Got "
                "(penalty={})".format(self.penalty)
            )
        if self.penalty == "none":
            if self.C != 1.0:  # default values
                warnings.warn(
                    "Setting penalty='none' will ignore the C and l1_ratio parameters"
                )
                # Note that check for l1_ratio is done right above
            C_ = np.inf
            penalty = "l2"
        else:
            C_ = self.C
            penalty = self.penalty
        if not isinstance(self.max_iter, _logistic.numbers.Number) or self.max_iter < 0:
            raise ValueError(
                "Maximum number of iteration must be positive; got (max_iter=%r)"
                % self.max_iter
            )
        if not isinstance(self.tol, _logistic.numbers.Number) or self.tol < 0:
            raise ValueError(
                "Tolerance for stopping criteria must be positive; got (tol=%r)"
                % self.tol
            )
        return C_, penalty   

    def set_classes(self_classes):
        n_classes = len(self_classes)
        classes_ = self_classes
        if n_classes < 2:
            raise ValueError(
                "This solver needs samples of at least 2 classes"
                " in the data, but the data contains only one"
                " class: %r"
                % classes_[0]
            )

        if len(self_classes) == 2:
            n_classes = 1
            classes_ = classes_[1:]
        
        return n_classes, classes_

    def set_fold_coef(self, X, y, C_, classes_, warm_start_coef, prefer, max_squared_sum, multi_class, solver, penalty, sample_weight, n_threads):
        path_func = _logistic.delayed(_logistic._logistic_regression_path)
        fold_coefs_ = _logistic.Parallel(n_jobs=self.n_jobs, verbose=self.verbose, prefer=prefer)(
            path_func(
                X,
                y,
                pos_class=class_,
                Cs=[C_],
                l1_ratio=self.l1_ratio,
                fit_intercept=self.fit_intercept,
                tol=self.tol,
                verbose=self.verbose,
                solver=solver,
                multi_class=multi_class,
                max_iter=self.max_iter,
                class_weight=self.class_weight,
                check_input=False,
                random_state=self.random_state,
                coef=warm_start_coef_,
                penalty=penalty,
                max_squared_sum=max_squared_sum,
                sample_weight=sample_weight,
                n_threads=n_threads,
            )
            for class_, warm_start_coef_ in zip(classes_, warm_start_coef)
        )

        return fold_coefs_

    def fit(self, X, y, sample_weight=None):
        solver = _logistic._check_solver(self.solver, self.penalty, self.dual)
        C_, penalty = self.warning_checks()

        if solver == "lbfgs":
            _dtype = np.float64
        else:
            _dtype = [np.float64, np.float32]

        X, y = self._validate_data(
            X,
            y,
            accept_sparse="csr",
            dtype=_dtype,
            order="C",
            accept_large_sparse=solver not in ["liblinear", "sag", "saga"],
        )
        _logistic.check_classification_targets(y)
        self.classes_ = np.unique(y)

        multi_class = _logistic._check_multi_class(self.multi_class, solver, len(self.classes_))

        # Good for small datasets - maybe don't focus as much on this
        if solver == "liblinear":
            if _logistic.effective_n_jobs(self.n_jobs) != 1:
                warnings.warn(
                    "'n_jobs' > 1 does not have any effect when"
                    " 'solver' is set to 'liblinear'. Got 'n_jobs'"
                    " = {}.".format(_logistic.effective_n_jobs(self.n_jobs))
                )
            self.coef_, self.intercept_, self.n_iter_ = _logistic._fit_liblinear(
                X,
                y,
                self.C,
                self.fit_intercept,
                self.intercept_scaling,
                self.class_weight,
                self.penalty,
                self.dual,
                self.verbose,
                self.max_iter,
                self.tol,
                self.random_state,
                sample_weight=sample_weight,
            )
            return self

        max_squared_sum = _logistic.row_norms(X, squared=True).max() if solver in ["sag", "saga"] else None
        n_classes, classes_ = self.set_classes(self.classes_)

        if self.warm_start:
            warm_start_coef = getattr(self, "coef_", None)
        else:
            warm_start_coef = None
        if warm_start_coef is not None and self.fit_intercept:
            warm_start_coef = np.append(
                warm_start_coef, self.intercept_[:, np.newaxis], axis=1
            )

        # Hack so that we iterate only once for the multinomial case.
        if multi_class == "multinomial":
            classes_ = [None]
            warm_start_coef = [warm_start_coef]
        if warm_start_coef is None:
            warm_start_coef = [None] * n_classes

        # The SAG solver releases the GIL so it's more efficient to use
        # threads for this solver.
        if solver in ["sag", "saga"]:
            prefer = "threads"
        else:
            prefer = "processes"

        n_threads = 1

        fold_coefs_ = self.set_fold_coef(X, y, C_, classes_, warm_start_coef, 
                                         prefer, max_squared_sum, multi_class, 
                                         solver, penalty, sample_weight, \
                                         n_threads)

        fold_coefs_, _, n_iter_ = zip(*fold_coefs_)
        self.n_iter_ = np.asarray(n_iter_, dtype=np.int32)[:, 0]

        n_features = X.shape[1]
        if multi_class == "multinomial":
            self.coef_ = fold_coefs_[0][0]
        else:
            self.coef_ = np.asarray(fold_coefs_)
            self.coef_ = self.coef_.reshape(
                n_classes, n_features + int(self.fit_intercept)
            )

        if self.fit_intercept:
            self.intercept_ = self.coef_[:, -1]
            self.coef_ = self.coef_[:, :-1]
        else:
            self.intercept_ = np.zeros(n_classes)

        return self