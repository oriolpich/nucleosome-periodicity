import pandas as pd


INCREASE_HEADER = ['nmuts_whole_nucleosome', 'observed_in',
              'observed_out', 'expected_in', 'expected_out', 'snr',
              'peak', 'cross_validation_max', 'empirical_pvalue_snr']

def load(file, field=None):
    df = pd.read_csv(file, sep='\t', names=INCREASE_HEADER)
    if field is None:
        return df
    else:
        return df[field]


def load_d(files):
    all = []
    for name, file in files.items():
        df = load(file)
        df['name'] = name
        all.append(df)
    return pd.concat(all, ignore_index=True)
