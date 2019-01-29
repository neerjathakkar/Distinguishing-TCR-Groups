# By Neerja Thakkar for "Balancing sensitivity and specificity in
# distinguishing TCR groups by CDR sequence similarity"
# See README for license information

# Smith Waterman scoring method

import swalign

# uses a gap penalty of 10
def getScoreSW(aa1, aa2, gap_penalty=-10):
    # set blosum scoring matrix
    scoring = swalign.ScoringMatrix('blosum_45.txt')

    # align the sequences
    sw = swalign.LocalAlignment(scoring, gap_penalty)
    alignment = sw.align(aa1, aa2)

    return alignment.score

# returns a normalized (between 0 and 1) SW distance
def getDistanceSW(aa1, aa2, length_dep=True, gap_penalty=-10):
    if length_dep:
        if len(aa1) != len(aa2):
            return 1.0
            
     # get self-scores for both
    self_score_1 = getScoreSW(aa1, aa1, gap_penalty)
    self_score_2 = getScoreSW(aa2, aa2, gap_penalty)

    # we want to divide by the min distance
    min_self_score = float(min(self_score_1, self_score_2))

    if min_self_score <= 0:
        print("warning: negative or zero self-score")
        return 1.0

    normalized = getScoreSW(aa1, aa2, gap_penalty)/min_self_score

    return 1.0 - normalized
