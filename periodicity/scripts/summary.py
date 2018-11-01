import click
import numpy as np
import pandas as pd
from scipy.stats import chi2, mannwhitneyu


def log_empirical_pvalue(val, sample):
    n = sum([1 for s in sample if s > val])
    try:
        p = n / len(sample)
        if p == 0:
            return 2
        else:
            return -np.log10(p)
    except ZeroDivisionError as e:
        return np.nan


def compute_summary(observed_chunks_file, random_chunks_file):
    """
    Generates summary spectral analysis statistics from observed and randomized 1Mb chunks
    """

    df_observed = pd.read_csv(observed_chunks_file, sep='\t')
    df_random = pd.read_csv(random_chunks_file, sep='\t')

    g = df_random.groupby(by='chunk_id')

    df_list = []
    for i in range(5, 20):
        df_list.append(pd.DataFrame(g['power_{0}p'.format(str(i))].median()))
        df_list.append(pd.DataFrame(g['power_{0}p'.format(str(i))].std()))
        df_list.append(pd.DataFrame(g['snr_{0}p'.format(str(i))].median()))
        df_list.append(pd.DataFrame(g['snr_{0}p'.format(str(i))].std()))

    table_fig7 = pd.concat(df_list, axis=1)
    cols = []
    for i in range(5, 20):
        cols += list(map(lambda x: x + str(i) + 'p', ['median_power_', 'std_power_', 'median_snr_', 'std_snr_']))
    table_fig7.columns = cols

    cols = ['chunk_id']

    for i in range(5, 20):
        cols += ['logpval_power_{0}p'.format(str(i)), 'logpval_snr_{0}p'.format(str(i))]

    pvals = pd.DataFrame(columns=cols)
    pvals_power = {str(i): [] for i in range(5, 20)}
    pvals_snr = {str(i): [] for i in range(5, 20)}

    chunks = []
    for i, chunk in enumerate(df_observed['chunk_id'].values):
        for j in range(5, 20):
            pvals_power[str(j)].append(log_empirical_pvalue(df_observed.loc[i, 'power_{0}p'.format(str(j))],
                                                            df_random[df_random['chunk_id'] == chunk][
                                                                'power_{0}p'.format(str(j))].values))

            pvals_snr[str(j)].append(log_empirical_pvalue(df_observed.loc[i, 'snr_{0}p'.format(str(j))],
                                                          df_random[df_random['chunk_id'] == chunk][
                                                              'snr_{0}p'.format(str(j))].values))

        chunks.append(chunk)

    pvals['chunk_id'] = np.array(chunks)
    for i in range(5, 20):
        pvals['logpval_power_{0}p'.format(str(i))] = np.array(pvals_power[str(i)])
        pvals['logpval_snr_{0}p'.format(str(i))] = np.array(pvals_snr[str(i)])

    table_fig7 = pd.merge(df_observed, table_fig7, right_on='chunk_id', left_on='chunk_id', right_index=True)
    table_fig7 = pd.merge(table_fig7, pvals, right_on='chunk_id', left_on='chunk_id')

    for i in range(5, 20):
        table_fig7['fold_power_increase_{0}p'.format(str(i))] = (table_fig7['power_{0}p'.format(str(i))] - table_fig7[
            'median_power_{0}p'.format(str(i))]) / table_fig7['median_power_{0}p'.format(str(i))]
        table_fig7['fold_snr_increase_{0}p'.format(str(i))] = (table_fig7['snr_{0}p'.format(str(i))] - table_fig7[
            'median_snr_{0}p'.format(str(i))]) / table_fig7['median_snr_{0}p'.format(str(i))]

    return table_fig7


def gtest(obs, exp):
    """

    Args:
        observed: array: observed values
        expected: array: expected values

    Returns:
        Returns the p-value obtained by the G-test

    """
    x = [val * np.log(val / exp[i]) for i, val in enumerate(obs)]
    chi = chi2(df=len(x) - 1)
    return chi.sf(2 * sum(x))


def componentwise_gtest(observed, expected):
    """

    Args:
        observed: array of float (observed freq for each bin)
        expected: array of float (expected freq for each bin)

    Returns:
        qvalues of every G-test consisting of one bin vs sum of the rest

    """

    n = len(observed)
    pvals = []
    log_odds_ratios = []
    obs_rest = sum(observed)
    exp_rest = sum(expected)
    for i in range(n):
        obs_rest -= observed[i]
        exp_rest -= expected[i]
        a = np.array([observed[i], obs_rest])
        b = np.array([expected[i], exp_rest])
        pvals.append(gtest(a, b))
        log_odds_ratios.append(np.log(observed[i] * exp_rest) - np.log(expected[i] * obs_rest))
        obs_rest += observed[i]
        exp_rest += expected[i]
    return pvals, log_odds_ratios


def odds_one_vs_rest(data):
    """
    data: one-dimensional array
    """

    data = np.array(data)
    expected = np.mean(data) * np.ones(len(data))
    nan_mask = np.isnan(data)
    data[nan_mask] = 1 / len(data)
    return componentwise_gtest(data, expected)


def get_metrics(table_fig7, reference_period):
    # index of reference period
    index_ref_period = list(range(6, 20)).index(reference_period)

    # snr in all chunks

    snrs = []
    df = table_fig7[table_fig7['maxp'] < 19.5]
    df = df[df['maxp'] > 5.5]
    for i in range(6, 20):
        snrs.append(df['snr_{0}p'.format(str(i))].values)

    x = np.concatenate([snrs[index_ref_period]], axis=0)
    y = np.concatenate(snrs[:index_ref_period] + snrs[index_ref_period + 1:], axis=0)
    _, pval_snr_all_chunks = mannwhitneyu(x, y)
    effect_snr_all_chunks = (np.median(x) - np.median(y)) / np.median(y)
    max_mean_snr_all_chunks = range(6, 20)[np.argmax(list(map(lambda x: 0 if len(x) == 0 else np.mean(x), snrs)))]

    # snr in chunks with peak period

    snrs = []
    for i in range(6, 20):
        df = table_fig7[table_fig7['logpval_snr_{0}p'.format(str(i))] == 2]
        df = df[df['maxp'] < i + 0.5]
        df = df[df['maxp'] > i - 0.5]
        snrs.append(df['snr_{0}p'.format(str(i))].values)
    x = np.concatenate([snrs[index_ref_period]], axis=0)
    y = np.concatenate(snrs[:index_ref_period] + snrs[index_ref_period + 1:], axis=0)
    _, pval_snr_peaks = mannwhitneyu(x, y)
    effect_snr_peaks = (np.median(x) - np.median(y)) / np.median(y)
    max_mean_snr_peaks = range(6, 20)[np.argmax(list(map(lambda x: 0 if len(x) == 0 else np.mean(x), snrs)))]

    # fold power increase in all chunks

    power = []
    df = table_fig7[table_fig7['maxp'] < 19.5]
    df = df[df['maxp'] > 5.5]
    for i in range(6, 20):
        power.append(df['fold_power_increase_{0}p'.format(str(i))].values)
    x = np.concatenate([power[index_ref_period]], axis=0)
    y = np.concatenate(power[:index_ref_period] + power[index_ref_period + 1:], axis=0)
    _, pval_fold_power_all_chunks = mannwhitneyu(x, y)
    effect_fold_power_all_chunks = (np.median(x) - np.median(y)) / np.median(y)
    max_mean_fold_power_all_chunks = range(6, 20)[np.argmax(list(map(lambda x: 0 if len(x) == 0 else np.mean(x), power)))]

    # fold power increase at peak period

    power = []
    for i in range(6, 20):
        df = table_fig7[table_fig7['maxp'] < i + 0.5].copy()
        df = df[df['maxp'] > i - 0.5]
        power.append(df['fold_power_increase_{0}p'.format(str(i))].values)
    x = np.concatenate([power[index_ref_period]], axis=0)
    y = np.concatenate(power[:index_ref_period] + power[index_ref_period + 1:], axis=0)
    _, pval_fold_power_peaks = mannwhitneyu(x, y)
    effect_fold_power_peaks = (np.median(x) - np.median(y)) / np.median(y)
    max_mean_fold_power_peaks = range(6, 20)[np.argmax(list(map(lambda x: 0 if len(x) == 0 else np.mean(x), power)))]

    # fold snr increase in all chunks

    snrs = []
    df = table_fig7[table_fig7['maxp'] < 19.5]
    df = df[df['maxp'] > 5.5]
    for i in range(6, 20):
        snrs.append(df['fold_snr_increase_{0}p'.format(str(i))].values)
    x = np.concatenate([snrs[index_ref_period]], axis=0)
    y = np.concatenate(snrs[:index_ref_period] + snrs[index_ref_period + 1:], axis=0)
    _, pval_fold_snr_all_chunks = mannwhitneyu(x, y)
    effect_fold_snr_all_chunks = (np.median(x) - np.median(y)) / np.median(y)
    max_mean_fold_snr_all_chunks = range(6, 20)[np.argmax(list(map(lambda x: 0 if len(x) == 0 else np.mean(x), snrs)))]

    # fold snr increase at peaks

    snrs = []
    for i in range(6, 20):
        df = table_fig7[table_fig7['maxp'] < i + 0.5]
        df = df[df['maxp'] > i - 0.5]
        snrs.append(df['fold_snr_increase_{0}p'.format(str(i))].values)
    x = np.concatenate([snrs[index_ref_period]], axis=0)
    y = np.concatenate(snrs[:index_ref_period] + snrs[index_ref_period + 1:], axis=0)
    ustat_fold_snr_peaks, pval_fold_snr_peaks = mannwhitneyu(x, y)
    effect_fold_snr_peaks = (np.median(x) - np.median(y)) / np.median(y)
    max_mean_fold_snr_peaks = range(6, 20)[np.argmax(list(map(lambda x: 0 if len(x) == 0 else np.mean(x), snrs)))]

    # significance power enrichment for each period

    discoveries = []
    for i in range(6, 20):
        discoveries.append(len(table_fig7[table_fig7['logpval_power_{0}p'.format(str(i))] == 2]))
    pvals_power_enrichment, logodds_power_enrichment = odds_one_vs_rest(discoveries)
    max_period_power_enrichment = range(6, 20)[np.argmax(logodds_power_enrichment)]

    # enrichemnt of periods in chunks significantly high in power

    discoveries = []
    for i in range(6, 20):
        discoveries.append(len(table_fig7[table_fig7['logpval_power_{0}p'.format(str(i))] == 2]))
    pval_power_enrichment, logodds_power_enrichment = odds_one_vs_rest(discoveries)
    max_logodds_power_enrichment = range(6, 20)[np.argmax(logodds_power_enrichment)]

    # enrichment of periods in chunks significantly high in snr

    discoveries = []
    for i in range(6, 20):
        discoveries.append(len(table_fig7[table_fig7['logpval_snr_{0}p'.format(str(i))] == 2]))
    pval_snr_enrichment, logodds_snr_enrichment = odds_one_vs_rest(discoveries)
    max_snr_enrichment = range(6, 20)[np.argmax(logodds_snr_enrichment)]

    # proportion of chunks at peak period

    df = table_fig7[table_fig7['maxp'] < reference_period + 0.5]
    df = df[df['maxp'] > reference_period - 0.5]
    proportion_max_ref_period = len(df) / len(table_fig7)

    return_tuple = (reference_period, effect_snr_all_chunks, pval_snr_all_chunks, max_mean_snr_all_chunks,
                    effect_snr_peaks, pval_snr_peaks, max_mean_snr_peaks,
                    effect_fold_power_all_chunks, pval_fold_power_all_chunks, max_mean_fold_power_all_chunks,
                    effect_fold_power_peaks, pval_fold_power_peaks, max_mean_fold_power_peaks,
                    effect_fold_snr_all_chunks, pval_fold_snr_all_chunks, max_mean_fold_snr_all_chunks,
                    effect_fold_snr_peaks, pval_fold_snr_peaks, max_mean_fold_snr_peaks,
                    pval_power_enrichment[index_ref_period], logodds_power_enrichment[index_ref_period],
                    max_logodds_power_enrichment,
                    pval_snr_enrichment[index_ref_period], logodds_snr_enrichment[index_ref_period], max_snr_enrichment,
                    proportion_max_ref_period)

    return return_tuple


def compute_metrics(df_summary):
    df_list = []
    cols = ['period',
            'effect_snr_all_chunks', 'pval_snr_all_chunks', 'max_mean_snr_all_chunks',
            'effect_snr_peaks', 'pval_snr_peaks', 'max_mean_snr_peaks',
            'effect_fold_power_all_chunks', 'pval_fold_power_all_chunks', 'max_mean_fold_power_all_chunks',
            'effect_fold_power_peaks', 'pval_fold_power_peaks', 'max_mean_fold_power_peaks',
            'effect_fold_snr_all_chunks', 'pval_fold_snr_all_chunks', 'max_mean_fold_snr_all_chunks',
            'effect_fold_snr_peaks', 'pval_fold_snr_peaks', 'max_mean_fold_snr_peaks',
            'pval_power_enrichment', 'logodds_power_enrichment', 'max_period_power_enrichment',
            'pval_snr_enrichment', 'logodds_snr_enrichment', 'max_period_snr_enrichment',
            'proportion_ref_period']

    for period in range(6, 20):
            d = dict(zip(cols, list(get_metrics(df_summary, period))))
            df_list.append(d)

    dh = pd.DataFrame(df_list)
    for col in cols:
        if col.startswith('pval'):
            dh['log_' + col] = dh[col].apply(lambda x: -np.log10(x))

    return dh


@click.command(context_settings={'help_option_names': ['h', '--help']})
@click.argument('observed', metavar='<OBS CHUNKS>', type=click.Path(exists=True))
@click.argument('random', metavar='<RAND CHUNKS>', type=click.Path(exists=True))
@click.argument('summary', metavar='<SUMMARY FILE>', type=click.Path())
@click.argument('metrics', metavar='<METRICS FILE>', type=click.Path())
def get_table(observed, random, summary, metrics):

    df_summary = compute_summary(observed, random)
    df_summary.to_csv(summary, sep='\t', index=False)

    df_metrics = compute_metrics(df_summary)
    df_metrics.to_csv(metrics, sep='\t', index=False)


if __name__ == '__main__':
    get_table()
