# -*- coding: utf-8 -*-

import csv
import heapq
import glob
import os.path
import sys
import traceback
from collections import defaultdict, namedtuple

import sn
import numpy as np


NaiveRecord = namedtuple("NaiveRecord", "ID REF ALT QUAL INFO")


def buggy(f):
    def inner(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except Exception as e:
            traceback.print_exc()
            import pdb; pdb.set_trace()

    return inner


required_fields = set(["AF", "ASN_AF", "AMR_AF", "AFR_AF", "EUR_AF"])

def naive_vcf_reader(handle):
    for line in handle:
        if line.startswith("#"):
            continue

        (_chrom, _pos, name,
         ref, alt, qual,
         _filter, raw_info) = line.split()[:8]

        info, non_snp = {}, False
        for kv in raw_info.split(";"):
            if "=" not in kv:
                continue

            k, v = kv.split("=")

            if k == "VT" and v != "SNP":
                non_snp = True

            if k in required_fields:
                info[k] = float(v)

        if non_snp:
            continue

        yield NaiveRecord(name, ref, alt.split(","), float(qual), info)


# Task #1
# Используя ... определите вероятность того, что вы американец (AMR),
# азиат (ASN), африканец (AFR), европеец (EUR). Должно получиться 4
# значения вероятностей.
def race(records):
    race_freqs = np.array([.25, .25, .25, .25], np.float64)
    races, eps = ["AMR_AF", "ASN_AF", "AFR_AF", "EUR_AF"], 10e-8
    for record in records.itervalues():
        snp_freq = np.array(
            [record.INFO.get(race_id, eps) for race_id in races],
            np.float64)

        race_freqs *= snp_freq * record.QUAL
        race_freqs /= sum(race_freqs)

    return race_freqs


# Task #2
# Используя данные ~152 генотипов из http://opensnp.org/dump_download
# определите 10 генетически ближайших к вам индивидуумов (например, по
# количеству совпадающих гаплотипов в снипах).
@buggy
def neighbours(snps, records, candidates, n=10):
    def rank(candidate):
        total = 0.

        try:
            for snp in sn.parse(candidate):
                for g in snp.genotype:
                    # Similar to 'lookup_snps' we treat alleles independantly.
                    if g in (snp.name, g) not in known:
                        continue

                    total += records[snp.name].INFO.get("AF", 0.)
        except Exception:
            return -1
        else:
            return total

    # *Both* alleles shoudl match!
    known = set((snp.id, snp.genotype) for snp in snps)
    return sorted(heapq.nlargest(n, candidates, rank))

# Task #3
# Используя 10 ближайших "родственников" из пункта 2, определите свой
# фенотип (например, взяв мажорирующие признаки из фенотипов
# "родственников").
@buggy
def phenotype(neighbours, traits):
    traits = dict((p["user_id"], p)
                  for p in csv.DictReader(open(traits), delimiter=";"))

    def inner():
        acc = defaultdict(lambda: defaultdict(int))

        for neighbour in map(os.path.basename, neighbours):
            user_id, _ = neighbour.split("_", 1)
            user_id = user_id[4:]  # len("user") == 4

            if user_id not in traits:
                continue

            for trait, option in traits[user_id].iteritems():
                if option == "-" or trait in ["user_id", "date_of_birth"]:
                    continue

        	acc[trait][option] += 1

        for trait, options in acc.iteritems():
            yield trait, max(options, key=lambda o: options[o])

    return sorted(inner())


def lookup_snps(snps, dbSNP):
    # Unfortuntately, dbSNP doesn't provide heterozygosity information,
    # (at least in the VCF format), so we treat two alleles
    # independantly.
    known = set((snp.name, g) for snp in snps for g in snp.genotype)
    resolved = {}

    for record in naive_vcf_reader(open(dbSNP)):
        if all((record.ID, str(alt)) not in known for alt in record.ALT):
            continue

        resolved[record.ID] = record

    return resolved


if __name__ == "__main__":
    try:
        dbSNP, opensnp, someone = sys.argv[1:]
    except ValueError:
        sys.exit("{0} dbSNP OPENSNP SOMEONE".format(__file__))

    print("Looking up sample SNPs in dbSNP ...")
    snps = sn.parse(someone)
    records = lookup_snps(snps, dbSNP)
    print("Done.")

    print("Task #1: race probabilities")
    print(race(records))

    print("Task #2: opensnp neighbours")
    found = neighbours(snps,
                       records,
                       glob.glob(os.path.join(opensnp, "user*.txt")))
    print(found)

    print("Task #3: phenotype")
    [phenotypes] = glob.glob(os.path.join(opensnp, "phenotypes_*.csv"))
    print(phenotype(found, phenotypes))
