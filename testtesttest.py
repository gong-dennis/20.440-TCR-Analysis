import pandas as pd

def main():
    for i in range(1,10):
        good_one = pd.read_csv(f'/Users/thomasusherwood/courses_code/20_440/20.440-TCR-Analysis/data/processed/cancer_annotated/LivMet_0{i}.csv')
        bad_one = pd.read_csv(f'/Users/thomasusherwood/courses_code/20_440/20.440-TCR-Analysis/data/processed/shared_seqs_deepcat/LivMet_0{i}.csv')
        good_one_set = set(good_one['aminoAcid'])
        bad_one_set = set(bad_one['aa_seqs'])
        good_one_set = bad_one_set
        good_one_set = set(good_one['aminoAcid'])
        print(good_one_set == bad_one_set)


if __name__ == '__main__':
    main()