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

