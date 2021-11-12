import pytest
from gd_ilp import run_optimizer

# some sample trees
stree = "(((a,b),c),(d,e))"
gt1 = "(((a,a),a),(a,a))"  # gene tree
gt2 = "(((b,b),c),d)"
gt3 = "(((a,a),a),(d,d))"
gt4 = "(((a,a),b),(a,d))"
gt5 = "(((b,(b,b)),c),d)"
gt6 = "(((a,a),a),a)"
gt7 = "((a,b),(a,d))"
gt8 = "(((a,e),(e,e)),((a,e),(a,a)))"
gt9 = "(((a,e),(e,e)),(c,(a,a)))"
gt10 = "(((b,(b,b)),((c,b),(c,b))),d)"
gt11 = "((b,(b,b)),(((c,b),(c,b)),d))"
gt12 = "((d,(b,b)),(((c,b),(c,b)),d))"
gt13 = "((((a,a),a),a),d)"
gfA = [gt6, gt2, gt3, gt4]  # 9; 4;
d = {0: [{0, 3}], 1: [], 2: [{2}], 3: [{0, 1, 2}], 4: [], 5: [], 6: [], 7: [{0, 2, 3}], 8: []}
# three comparable duplications from gt6, are nicely distributed to fixed locations
gfB = [gt1, gt2, gt3, gt4]  # 10; 4;     1, 0, 0, 1, 0, 1, 0, 1, 0
# # gt1 checks if proper chain is considered instead of sum of comparable duplications
gfC = [gt1, gt5, gt3, gt7]  # 10; 4;     1, 0, 0, 2, 0, 1, 0, 0, 0
d1 = {0: [{0, 3}], 1: [], 2: [], 3: [{0, 1, 2}, {0, 1, 2}], 4: [], 5: [{2}], 6: [], 7: [], 8: []}
# # checks level 2 chain in optimal solution;
# # gt7 replaces gt4 so there is no longer forced duplication at leaf a;
# # so chain from gt1, should mix with now longer chain from gt5 (gt5 instead of gt2)
gfD = [gt8, gt9]
gfE = [gt10, gt11, gt12]
gfF1 = [gt11, gt13]
gfF2 = [gt10, gt13]
# # (((a7,b8)3,c4)1,(d5,e6)2)0 ;
# # F1 - st node 1 has chain [0,1,2] but duplication gt node 0 is not mapped to this node in optimal solution
# # correct; other constraints guarantee that 0 cannot be mapped below root
gf = [gt1, gt5, gt3, gt4]  # 11; 5;     1, 0, 0, 2, 0, 1, 0, 1, 0
gt14 = "((a,a),b)"
gt15 = "((c,c),c)"
gfG = [gt14, gt15]
gt16 = "((((((a,a),a),b),(a,b)),c),(c,(c,(c,(c,(c,c))))))"
gt17 = "((((((a,a),a),b),(a,b)),c),(c,(c,(c,(c,c)))))"

gt18 = "(((a,a),a),((b,b),b))"
gt19 = "((a,b),c)"
st1 = "(a,(b,c))"
st2 = "(b,(a,c))"

gt20 = "((((SCHPO,SCHPO),YEAST),((((((((((((HUMAN,PANTR),MACMU),(BOVIN,CANFA)),(MOUSE,RAT)),MONDO),CHICK),XENTR),(BRARE,(TETNG,GASAC))),CIOIN),((DROME,DROME),(AEDAE,ANOGA))),SCHMA),(CAEEL,CAEBR))),ARATH)"
st_treefam = "((((((((AEDAE,ANOGA),(DROME,DROPS)),APIME),SCHMA),(((((((BOVIN,CANFA),(((HUMAN,PANTR),MACMU),(MOUSE,RAT))),MONDO),CHICK),XENTR),(BRARE,((TETNG,FUGRU),(ORYLA,GASAC)))),CIOIN)),(CAEBR,CAEEL)),(YEAST,SCHPO)),(ARATH,ORYSA))"
d2 = {0: [], 1: [], 2: [], 3: [], 4: [], 5: [], 6: [], 7: [{0}], 8: [], 9: [], 10: [{0}], 11: [], 12: [], 13: [], 14: [], 15: [], 16: [], 17: [], 18: [], 19: [], 20: [], 21: [], 22: [], 23: [], 24: [], 25: [], 26: [], 27: [], 28: [], 29: [], 30: [], 31: [{0}], 32: [], 33: [], 34: [], 35: [], 36: [], 37: [{0}], 38: [], 39: [], 40: [], 41: [], 42: [], 43: [], 44: [], 45: [], 46: [], 47: [], 48: [], 49: [], 50: [], 51: [], 52: [], 53: [], 54: []}


@pytest.mark.parametrize("params,totalscore,gtsupport", [
    ([stree, [gt1]], 3.0, {}),               # test computation of ME score - longest chain = 3
    ([stree, [gt1, gt2]], 3.0, {}),          # test movement to node (a,b)
    ([stree, [gt1, gt2], 1], 4.0, {}),       # interval_limit 1 prevents from movement of duplications
    ([stree, gfA], 3.0, {}),                 # unrestricted model all mapped to root
    ([stree, gfA, 0, 0], 4.0, {}),           # max_converted_spec set to 0, restricted model, top duplication from gt6
    #                                          at the root; 1, 0, 0, 1, 0, 1, 0, 1, 0
    ([stree, gfA, 0, 0], 4.0, d),            # to test support for duplication episodes
    # top duplication from gt6, and gt4 are mapped in solution to root S; hence d[0] = [{0,3}]
    ([stree, gfB, 0, 0], 4.0, {}),
    ([stree, gfC, 0, 0], 4.0, d1),           # tests if partial score = max chain length instead of sum of all
    # duplications; also checks support for a chain
    ([stree, gfD, 0, 0], 4.0, {}),           # check equal maximal chains
    ([stree, gfE, 0, 0], 4.0, {}),           # check multiple chains not in the root, intervals with common one node
    ([stree, gfF1, 0, 0], 4.0, {}),
    ([stree, gfF2, 0, 0], 3.0, {}),
    ([stree, gf, 0, 0], 5.0, {}),
    ([stree, gfG], 2.0, {}),
    ([stree, gfG, 0, 0], 3.0, {}),
    ([stree, [gt16]], 6.0, {}),              # chain of speciations and duplications, two speciations -> duplications
    # should be changed for optimal solution
    ([stree, [gt16], 0, 0], 9.0, {}),
    ([stree, [gt17]], 6.0, {}),
    ([stree, [gt17], 0, 0], 8.0, {}),
    ([st1, [gt18]], 3.0, {}),
    ([st2, [gt18]], 3.0, {}),
    ([gt19, [gt18]], 3.0, {}),
    ([st1, [gt18, gt19]], 3.0, {}),
    ([st2, [gt18, gt19]], 3.0, {}),
    ([st1, [st1, st1]], 0.0, {}),
    ([st1, [st1, gt18]], 3.0, {}),
    ([st_treefam, [gt20],1,0], 4.0, d2)  # TreeFam tree 31, lca-scenario
])


def test_run_optimizer(params, totalscore, gtsupport):
    (score, support) = run_optimizer(*params)
    assert score == totalscore
    if not gtsupport == {}:                  # test support only for some tests
        assert support == gtsupport

#  --time_limit 40 tests/treefam.txt => 25
#  --time_limit 20 tests/treefam.txt => no solution => 0
