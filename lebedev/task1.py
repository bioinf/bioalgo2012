# -*- coding: utf-8 -*-

# from __future__ import print_function

import csv
import itertools
import errno
import os
import os.path
import sys
import tempfile
import time
from collections import defaultdict

import numpy as np
from Bio import SeqIO
from Bio.Data import IUPACData
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from sklearn import ensemble
from sklearn.cross_validation import train_test_split
from sklearn.externals import joblib
from sklearn.utils import fixes


def GC_skew(seq, window=100):
    """Calculates GC skew (G-C)/(G+C) for multuple windows along the sequence.

    Returns a list of ratios (floats), controlled by the length of the sequence
    and the size of the window.

    Does NOT look at any ambiguous nucleotides.
    """
    values = []
    for i in range(0, len(seq), window):
        s = seq[i: i + window]
        g = s.count('G') + s.count('g')
        c = s.count('C') + s.count('c')

        # HACK(Sergei): default implementation raises 'ZeroDivisionError''
        # on non-GC sequences.
        if not g + c:
            values.append(0.)
        else:
            skew = (g - c) / float(g + c)
            values.append(skew)
    return values or [0.]


def molecular_weight(seq):
    """Calculate the molecular weight of a DNA sequence."""
    weight_table = IUPACData.unambiguous_dna_weights
    # HACK(Sergei): once again, fails with KeyError on 'N's.
    return np.sum(weight_table.get(x, 0.0) for x in seq)


chr_map = map(chr, xrange(256))
chr_map[ord("A")] = "0"
chr_map[ord("C")] = "1"
chr_map[ord("G")] = "2"
chr_map[ord("T")] = "3"
chr_map = "".join(chr_map)

def kmers(seq, k):
    seen = defaultdict(int)
    seq = str(seq)
    for i in xrange(0, len(seq) - k):
        chunk = seq[i:i + k]
        if "N" not in chunk:
            seen[int(chunk.translate(chr_map), base=4)] += 1

    q25 = np.percentile(seen.values(), .25)
    q75 = np.percentile(seen.values(), .75)
    result = []

    for kmer, count in seen.iteritems():
        if count < q25 or count > q75:
            continue

        result.append(kmer * count)
    else:
        return np.array(result)


# 4-element bins
bins = sorted(set("".join(bin) for bin in
                  itertools.product("ACGT", repeat=2)))

def hashed_kmers(seq, k=2):
    seen = dict.fromkeys(bins, 0)
    seq = str(seq)
    for i in xrange(0, len(seq) - k):
        seq_chunk = seq[i:i + k]
        if "N" not in seq_chunk:
            seen[seq_chunk] += 1

    kmers = np.zeros(len(bins))
    for idx, count in enumerate(seen.itervalues()):
        kmers[idx] = count

    return kmers


memory = joblib.Memory(cachedir=tempfile.gettempdir(), verbose=0)

@memory.cache
def generate_interesting_repeats(k):
    result = []

    if k % 2:
        chunks = generate_interesting_repeats(k // 2)
    else:
        chunks = []

    for n in ["A", "C", "G", "T"]:
        result.append("".join(itertools.repeat(n, k)))  # ex.: AAAAAA

        for chunk in chunks:
            result.append(chunk + n + chunk)            # ex.: AAACAAA

    return result


def nucleotide_repeats(seq, k=8):
    return np.fromiter(
        (seq.count(repeat) for repeat in generate_interesting_repeats(k)),
        np.int64)


class Genome(object):
    cache = {}

    def __init__(self, base_dir):
        self.base_dir = base_dir

    def features(self, chrom, start, end):
        if chrom not in self.cache:
            path = os.path.join(self.base_dir, chrom + ".fasta")
            self.cache[chrom], = SeqIO.parse(open(path), "fasta")

        chrom = self.cache.get(chrom).seq
        seq = chrom[start:end]
        protein = seq.translate()

        # XXX this is obviously BS, since most of the sequences don't look
        # like proper peptides, but hey, why not?
        analysis = ProteinAnalysis(str(protein))

        helix, turn, sheet = analysis.secondary_structure_fraction()
        trivial = [start, end, end - start,
                   hash(str(chrom)),
                   molecular_weight(seq),
                   protein.count("M"), protein.count("*"),
                   analysis.aromaticity(),
                   helix, turn, sheet,
                   np.median(GC_skew(seq))]

        return np.concatenate([trivial,
                               np.array(analysis.get_amino_acids_percent().values()),
                               hashed_kmers(seq),
                               nucleotide_repeats(seq, k=5)])


def extract_features(genome, path):
    print("Extracting features from {0} ...".format(path))

    chrom = os.path.basename(path)
    data  = (map(int, r) for r in csv.reader(open(path), delimiter=" "))
    features, labels = [], []
    for row in data:
        start, end = row[:2]
        features.append(genome.features(chrom, start, end))

        if len(row) == 3:
            labels.append(row[-1])

    return np.array(features), np.array(labels)


def cutoff_bootstrap(features, labels, n=50):
    freqs = fixes.Counter(labels)
    seen  = defaultdict(int)

    cutoff = np.percentile(freqs.values(), 0.25)
    print("Estimated cutoff: {0}.".format(cutoff))

    new_features, new_labels = [], []
    for f, l in itertools.izip(features, labels):
        if freqs[l] < cutoff:
            print("Dropping class {0}, freq = {1}.".format(l, freqs[l]))
            continue
        elif seen[l] == n:
            continue

        seen[l] += 1
        new_features.append(f)
        new_labels.append(l)

    return np.array(new_features), np.array(new_labels)


def get_classifier(genome, chroms, force=False, n_jobs=4):
    pickled_clf = ".classifier"
    pickled_features = "all_features.dump"
    if not force and os.path.exists(pickled_clf):
        print("Loading classifier from {0} ...".format(pickled_clf))
        clf = joblib.load(pickled_clf)
    else:
        # ~43.3900%
        clf = ensemble.RandomForestClassifier(
            min_samples_leaf=15, min_samples_split=150, compute_importances=True,
            n_jobs=4)

        if os.path.exists(pickled_features):
            print("Loading features from {0} ...".format(pickled_features))
            features, labels = joblib.load(pickled_features)
        else:
            features, labels = extract_all_features(genome, chroms, n_jobs=n_jobs)
            joblib.dump((features, labels), pickled_features)

        features_tr, features_te, labels_tr, labels_te = \
            train_test_split(features, labels, test_size=0.4)

        print("Training {0} ...".format(clf))
        clf.fit(features, labels)
        print("Done, saving params ...")
        joblib.dump(clf, pickled_clf, compress=9)

        print("Evaluating fit ...")
        print("Test score: ", clf.score(features_te, labels_te))

    return clf


def extract_all_features(genome, chroms, n_jobs=4):
    chunks = joblib.Parallel(n_jobs=n_jobs)(
        joblib.delayed(extract_features)(genome, chrom) for chrom in chroms)

    all_features, all_labels = [], []
    for features, labels in chunks:
        all_features.append(features)
        all_labels.append(labels)

    features = np.concatenate(all_features)
    labels   = np.concatenate(all_labels)

    print("Cleaning up training data ...")
    features_tr, labels_tr = cutoff_bootstrap(features, labels, n=100)

    print("Done.")
    return features, labels


TIMESTAMP = time.time()

def predict_all(clf, chroms, n_jobs=4):
    chunks = joblib.Parallel(n_jobs=n_jobs)(
        joblib.delayed(extract_features)(genome,
            os.path.join(predict_dir, chrom)) for chrom in chroms)

    for chrom, (features, _labels) in itertools.izip(chroms, chunks):
        predict(clf, chrom, features)


def predict(clf, chrom, features):
    result_dir = "result-{0:.8n}".format(TIMESTAMP)

    try:
        os.mkdir(result_dir, 0o755)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    positions = []
    for f in features:
        # HACK(Sergei): oopsy, abstraction leak.
        start, end = map(int, f[:2])
        positions.append((start, end))

    print("Predicting {0} ...".format(chrom))
    labels = clf.predict(features)

    with open(os.path.join(result_dir, chrom), "w") as result:
        for (start, end), label in itertools.izip(positions, labels):
            result.write("{0} {1} {2}\n".format(start, end, label))

        result.flush()


if __name__ == "__main__":
    try:
        root_dir, = sys.argv[1:]
    except ValueError:
        print("usage: {0} path/to/data/root".format(__file__))
        sys.exit(1)
    else:
        # XXX the expected root layout is:
        # root/genome/chrX,
        # root/train/chrX,
        # root/test/chrX,
        # where 'X' is a chromosome number -- 1-22, X, Y.
        genome_dir, train_dir, predict_dir =  \
            [os.path.join(root_dir, f) for f in ["genome", "train", "test"]]

    genome = Genome(genome_dir)
    chroms = os.listdir(train_dir)

    clf = get_classifier(genome,
        [os.path.join(train_dir, chrom) for chrom in chroms],
        n_jobs=4)

    predict_all(clf, chroms, n_jobs=8)
