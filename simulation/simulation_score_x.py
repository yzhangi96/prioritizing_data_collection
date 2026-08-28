import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def compute_score_x(reference, source, random_state=None):
    reference = np.asarray(reference, dtype=float)
    source = np.asarray(source, dtype=float)
    train_n = min(len(reference), len(source) // 2)
    if train_n < 1:
        raise ValueError("Score X requires training and validation rows.")

    rng = np.random.RandomState(random_state)
    reference_index = rng.choice(len(reference), size=train_n, replace=False)
    source_index = rng.permutation(len(source))
    source_train = source[source_index[:train_n]]
    source_validation = source[source_index[train_n:]]

    classifier = make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=10000, tol=1e-4),
    )
    classifier.fit(
        np.concatenate([source_train, reference[reference_index]]),
        np.concatenate([np.ones(train_n), np.zeros(train_n)]),
    )
    return float(classifier.predict_proba(source_validation)[:, 1].mean())
