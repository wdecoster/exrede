import pandas as pd
from argparse import ArgumentParser
from scipy import stats


def main():
    args = get_args()
    df = pd.read_csv(args.insertion_excess, sep="\t", header=None, names=['locus', 'excess'])
    print(
        df[(stats.zscore(df.excess) > 3)].to_csv(sep="\t", index=False),
        end="")


def get_args():
    parser = ArgumentParser(description="Find outliers in insertion excess "
                            "as quantified by windowed-insertion-excess.py")
    parser.add_argument("insertion_excess")
    return parser.parse_args()


if __name__ == '__main__':
    main()
