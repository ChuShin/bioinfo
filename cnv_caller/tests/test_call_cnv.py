import pytest
from cnv_caller.src.call_cnv import (load_covs,
                                     normalize,
                                     normalize_sample,
                                     call_he_events,
                                     print_he_output,
                                     find_seeds,
                                     extend_seeds,
                                     is_hit,
                                     call_he_type)


@pytest.mark.parametrize(
    'input,expected',
    [([], []),
     (['norm/norm', 'norm/del', 'norm/del', 'norm/del', 'norm/norm'], ['', 'norm/del', 'norm/del', 'norm/del', '']),
     (['norm/norm', '', 'norm/del', 'norm/del', 'norm/norm'], ['', '', '', '', '']),
     (['del/del', 'del/del', 'norm/del', 'norm/del', 'norm/norm'], ['', '', '', '', '']),
     ]
)
def test_find_seeds(input, expected):
    assert find_seeds(input) == expected

@pytest.mark.parametrize(
    'scoreA, scoreB, expected',
    [
        (-1, -1, 'norm/norm'), (0, 0, 'norm/norm'), (1, 1, 'norm/norm'),
        (-1.01, -1.01, 'del/del'), (-1.01, 1, 'del/norm'), (-1.01, 1.01, 'del/dup'),
        (-1, -1.01, 'norm/del'), (0, 1, 'norm/norm'), (1, 1.01, 'norm/dup'),
        (1.01, -1.01, 'dup/del'), (1.01, 0, 'dup/norm'), (1.01, 1.01, 'dup/dup'),
     ]
)
def test_call_he_type(scoreA, scoreB, expected):
    assert call_he_type(scoreA,scoreB) == expected

