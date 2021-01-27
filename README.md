# hicpileup

hicpileup generates HIC pileups around a BED file features

## Dependencies



## Usage

```
$ ./hicpileup.py
usage: hicpileup.py [-h] -c COOLER -cs CSIZES -fs [FEATURES [FEATURES ...]]
                    -fn [FEATURENAMES [FEATURENAMES ...]] -of OUTFILE
                    [-t THREADS]

Create HIC pileups around regions to aggregate

optional arguments:
  -h, --help            show this help message and exit
  -c COOLER, --cooler COOLER
                        coolerData
  -cs CSIZES, --csizes CSIZES
                        chrom-sizes file
  -fs [FEATURES [FEATURES ...]], --features [FEATURES [FEATURES ...]]
                        feature files
  -fn [FEATURENAMES [FEATURENAMES ...]], --featureNames [FEATURENAMES [FEATURENAMES ...]]
                        names for each feature file, used as plot labels
  -of OUTFILE, --outfile OUTFILE
                        output figure name
  -t THREADS, --threads THREADS
                        threads for processing pool
```


# Example output using two feature files:

![example](data/test1.png)


