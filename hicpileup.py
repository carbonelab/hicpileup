#!/usr/bin/env python3

"""
Creates HIC pileup plots using hic data, and a list of feature files.
"""

import bioframe
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pandas as pd
import cooler
import cooltools
import cooltools.expected
import cooltools.snipping
import os
import multiprocess
import numpy as np
import argparse
import warnings


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(prog="hicpileup.py",
                                     description="Create HIC pileups\
                                          around regions to aggregate")
    parser.add_argument("-c", "--cooler", help="coolerData", required=True)
    parser.add_argument("-cs", "--csizes", help="chrom-sizes file",
                        required=True)
    parser.add_argument("-fs", "--features",
                        help="feature files", nargs="*", required=True)
    parser.add_argument("-fn", "--featureNames", nargs="*", 
                        help="names for each feature file, used as plot labels", required=True)
    parser.add_argument("-of", "--outfile", help="output figure name",
                        required=True)
    parser.add_argument("-bs", "--binsize", help="cooler binsize",
                        required=False, default=10000, type=int)
    parser.add_argument("-t", "--threads", help="threads for processing pool",
                        required=False, default=5, type=int)
    return parser.parse_args()


def make_expected(clr, chromsizes, threads):
    cs = pd.read_csv(chromsizes, sep="\t", names=["chrom", "end"])
    cs['start'] = 0
    cs = cs[['chrom', 'start', 'end']]
    sprts = list(cs.itertuples(index=False, name=None))
    cs=bioframe.parse_regions(sprts)

    with multiprocess.Pool(threads) as pool:
        expected_df = cooltools.expected.diagsum(
            clr,
            regions=cs,
            transforms={
                'balanced': lambda p: p['count'] * p['weight1'] * p['weight2']
            },
            map=pool.map)

    expected_df['balanced.avg'] = \
        expected_df['balanced.sum'] / expected_df['n_valid']
    return expected_df


def snip_bed_features(bedfile, cs, binsize, window_halfsize):
    # organize features
    feat = pd.read_csv(bedfile,
                       sep="\t",
                       usecols=[0, 1, 2],
                       names=["chrom", "start", "end"])
    feat['mid'] = (feat['start'] + feat['end']) / 2

    # snip feature bins
    snipping_windows = cooltools.snipping.make_bin_aligned_windows(
        binsize,
        feat['chrom'],
        feat['mid'],
        window_halfsize)
    snipping_windows = cooltools.snipping.assign_regions(
        snipping_windows,
        cs)

    return snipping_windows


def make_pileup(csizes, expected_df, features, clr, threads):

    cs = pd.read_csv(csizes, sep="\t", names=["chrom", "end"])
    cs['start'] = 0
    cs = cs[['chrom', 'start', 'end']]

    # as lsit of tuples
    sprts = list(cs.itertuples(index=False, name=None))
    cs=bioframe.parse_regions(sprts)

    snippings = {}
    piles = {}
    for f in features:
        snippings[f] = snip_bed_features(features[f], cs, 10000, 1000000)
        oe_snipper = cooltools.snipping.ObsExpSnipper(clr,
                                                      expected=expected_df,
                                                      regions=cs)
        with multiprocess.Pool(threads) as pool:
            oe_pile = cooltools.snipping.pileup(
                snippings[f],
                oe_snipper.select, oe_snipper.snip,
                map=pool.map)
            piles[f] = np.nanmean(oe_pile, axis=2)
    return piles


def main():
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    args = parse_args()
    clr = args.cooler
    csizes = args.csizes
    features = args.features
    names = args.featureNames
    outfile = args.outfile
    threads = args.threads
    nconds = len(features)
    binsize = args.binsize
    conditions = names
    if len(conditions) != len(features):
        raise ValueError("The number of feature names does not equal the number of feature files, pleas specify one featureName per file") 
    clr = cooler.Cooler(clr)
    expected_df = make_expected(clr, csizes, threads)
    features = {os.path.basename(f): f for f in features}
    piles = make_pileup(csizes, expected_df, features, clr, threads)

    flank = 1000000
    gs = GridSpec(nrows=1, ncols=len(conditions) +
                  1, width_ratios=[20] * len(conditions) + [2]) 
    gs.update(wspace=0.25)
    plt.figure(figsize=(5 * len(conditions) + 2, 5))

    opts = dict(
        vmin=-0.5,
        vmax=0.5,
        extent=[-flank//1000, flank//1000, -flank//1000, flank//1000],
        cmap='coolwarm'
    )

    for i, cond in enumerate(features):
        ax = plt.subplot(gs[i])
        img = ax.matshow(
            np.log2(piles[cond]),
            **opts)
        ax.xaxis.tick_bottom()
        if i > 0:
            ax.yaxis.set_visible(False)
        plt.title(conditions[i])
        plt.xlabel('relative position, kbp')
        plt.ylabel('relative position, kbp')

    ax = plt.subplot(gs[nconds])
    plt.colorbar(img, cax=ax, label='log2 mean obs/exp')
    plt.suptitle(f"HIC pileup")
    plt.savefig(outfile)


if __name__ == "__main__":
    main()
