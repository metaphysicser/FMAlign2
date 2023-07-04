"""

Author:zhai1xiao
Create Date:2021-12-06
Update Date:2021-12-19

"""
import re
import argparse


def read_fasta(file_path):
    """
    读取fasta文件
    :param file_path: 路径
    :return: identifiers [] 序列id的列表
             sequences [] 序列的列表
    """
    identifiers = []
    sequences = []
    with open(file_path, 'r') as fin:
        curr = ""
        flag = False
        for line in fin:
            if line[0] == '>':
                if flag:
                    sequences.append(curr.upper())
                curr = ""
                flag = True
                identifiers.append(line[1:].strip())
            else:
                curr += line.strip()
        sequences.append(curr.upper())
    return identifiers, sequences


def preprocess(sequences):
    """
    处理fasta序列中的非法字符（非字母用N代替）

    :param sequences 待处理的序列列表
    :return: processed [] 处理后的序列列表
    """
    processed = []
    for i in sequences:
        tmp = re.sub("[^ACTGNactgn-]", "N", str(i).replace('U', 'T'))
        processed.append(tmp)
    return processed

def score_of(curr_clm: dict, matchS, mismatchS, gap1S, gap2S):
    """
    序列单列打分

    Args:
        sequences  待打分的序列
        matchS     match(碱基相同且非N)的得分
        mismatchS  mismatch(碱基非空非N，但不相同)的得分
        gap1S      gap比对碱基，N比对N，N比对碱基的得分
        gap2S      gap比对gap，gap比对N的得分

    Return:
        该列的sp score
    """
    match = 0  # nongap==nongap
    mismatch = 0  # nongap!=nongap
    gap1 = 0  # gap-nongap
    gap2 = 0
    A, C, G, T, N, dash = curr_clm['A'], curr_clm['C'], curr_clm['G'], curr_clm['T'], curr_clm['N'], curr_clm['-']
    match = (A * (A - 1) + C * (C - 1) + G * (G - 1) + T * (T - 1)) // 2
    mismatch = ((A + C) * (G + T)) + A * C + G * T
    gap1 = (A + C + G + T + N) * dash
    gap2 = (dash * (dash - 1) // 2) + ((N * (N - 1)) // 2) + (A + C + G + T) * N
    score = (match * matchS) + (mismatch * mismatchS) + (gap1 * gap1S) + (gap2 * gap2S)
    return score


def evaluate(sequences, matchS, mismatchS, gap1S, gap2S):
    """
    序列打分

    Args:
        sequences  待打分的序列
        matchS     match(碱基相同)的得分
        mismatchS  mismatch(碱基非空，但不相同)的得分
        gap1S      gap比对碱基的得分
        gap2S      gap比对gap的得分

    Return:
        sp score
    """
    print("calculating", end="", flush=True)
    score_list = []
    for j in range(len(sequences[0])):
        curr_clm: dict = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0, '-': 0}
        for i in range(len(sequences)):
            curr_clm[sequences[i][j]] += 1
        score_list.append(score_of(curr_clm, matchS, mismatchS, gap1S, gap2S))
        if (j + 1) % (len(sequences[0]) // 3) == 0: print('.', end="", flush=True)
    print("done")
    return sum(score_list)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, required=False, default='output.fmaligned2.fasta', help='Fasta file path to be scored.')
    parser.add_argument('--match', type=float, default=0, help='match score, default=1.')
    parser.add_argument('--mismatch', type=float, default=1, help='mismatch score, default=-1.')
    parser.add_argument('--gap1', type=float, default=2,help='gap-base score, default=-2.')
    parser.add_argument('--gap2', type=float, default=2, help='gap-gap score, default=0.')
    args = parser.parse_args()

    filepath = args.input
    matchScore = args.match
    mismatchScore = args.mismatch
    gap1Score = args.gap1
    gap2Score = args.gap2

    _, sequences = read_fasta(filepath)
    print(str(len(sequences)) + " sequences found")
    if len(sequences) <= 1: exit()

    processedSeq = preprocess(sequences)

    sp = evaluate(processedSeq, matchScore, mismatchScore, gap1Score, gap2Score)  # pure sp
    avgSp = sp / (len(sequences) * (len(sequences) - 1) // 2)  # avg sp
    scaledSP = avgSp / len(sequences[0])  # scaled sp
    print("SP score: " + str(sp))
    print("Avg SP score: " + str(avgSp))
    print("Scaled SP score: " + str(scaledSP))
    print("Finish!")