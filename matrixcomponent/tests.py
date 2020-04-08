import numpy as np
import pandas as pd

from matrixcomponent.utils import (
    find_groups,
    path_boundaries,
    path_dividers,
    self_loops,
    sort_and_drop_duplicates,
)


def test_path_boundaries():
    links = np.array([[0, 1], [1, 2], [2, 0]])
    assert np.array_equal(path_boundaries(links), [True, False, True])


def test_self_loops():
    links = np.array([[0, 1], [1, 1], [1, 2], [2, 1],
                      [1, 1], [1, 3], [3, 3], [3, 0]])
    loops = np.unique(links[self_loops(links)], axis=0)
    assert np.array_equal(loops, [[1, 1], [3, 3]])


def test_path_dividers():
    bins = [1, 3, 4, 8, 15]

    # test the simple case separately
    links = np.array([[3, 1], [4, 1]])
    assert np.array_equal(links[path_dividers(links, bins)], links)

    # test the general case
    links = np.array([
        [3, 1],  # u > d => yes
        [1, 3],  # no bins within [1 + 1, 3) range => no
        [3, 5],  # bin 4 is within [3 + 1, 5) range => yes
        [5, 9],  # bin 8 is within [5 + 1, 9) => yes
        [5, 8],  # no bins within [5 + 1, 9) => no
        [3, 4],  # empty bin range [3 + 1, 4) => no
        [4, 1],  # u > d => yes
    ])
    assert np.array_equal(links[path_dividers(links, bins)],
                          [[3, 1], [3, 5], [5, 9], [4, 1]])


def test_find_groups():
    data = np.array([
        [1, 2], [1, 2], [1, 3],
        [2, 1], [2, 1], [2, 1], [2, 2],
        [3, 3], [3, 3], [3, 4], [3, 4], [3, 5]
    ])
    assert np.array_equal(find_groups(data),
                          [(0, 2), (2, 3),
                           (3, 6), (6, 7),
                           (7, 9), (9, 11), (11, 12)])

    assert np.array_equal(find_groups(np.array([[]])), [])
    assert np.array_equal(find_groups(np.array([[1, 2]])), [(0, 1)])


def test_sort_and_drop_duplicates():
    df = pd.DataFrame.from_dict({
        "from":       [1, 3, 2, 3, 2, 0, 5, 4, 1, 2],
        "to":         [2, 2, 4, 2, 3, 1, 4, 3, 2, 3],
        "path_index": [0, 3, 1, 2, 2, 3, 2, 1, 3, 2],
    })

    expected = pd.DataFrame.from_dict({
        "from":       [0, 1, 1, 2, 2, 3, 3, 4, 5],
        "to":         [1, 2, 2, 3, 4, 2, 2, 3, 4],
        "path_index": [3, 0, 3, 2, 1, 2, 3, 1, 2],
    })  # only one duplicate (2, 3, 2)

    assert np.array_equal(sort_and_drop_duplicates(df).to_numpy(),
                          expected.to_numpy())