from ....algorythm.ga.selection import tournament_selection


def test_tournament_selection():
    p = [1, 0]
    assert tournament_selection(p) == p[1]
    p.append(-1)
    assert tournament_selection(p, selection_force=3) == p[-1]
