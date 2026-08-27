#!/usr/bin/env python3

import argparse
import gzip
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def open_file(filename):
    """Open normal or gzipped files."""
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    else:
        return open(filename)



def parse_ann(info):
    """
    Extract SnpEff ANN effects.
    Returns a list of effects.
    """
    effects = []
    for item in info.split(";"):
        if item.startswith("ANN="):
            annotations = item[4:].split(",")
            for ann in annotations:
                fields = ann.split("|")
                if len(fields) > 1:
                    effects.append(fields[1])
    return effects

def classify_effects(effects):
    # gene-associated effects
    gene_effects = {
        "synonymous_variant",
        "missense_variant",
        "stop_gained",
        "stop_lost",
        "start_lost",
        "frameshift_variant",
        "inframe_insertion",
        "inframe_deletion",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "splice_region_variant"
    }
    if len(effects) == 0:
        return "non_gene"
    # synonymous takes priority
    if "synonymous_variant" in effects:
        return "gene_synonymous"
    # any coding damaging effect
    if any(
        e in gene_effects
        for e in effects
    ):
        return "gene_nonsynonymous"
    return "non_gene"



def main():
    parser = argparse.ArgumentParser(
        description="Plot SnpEff variant density"
    )
    parser.add_argument(
        "-v",
        "--vcf",
        required=True,
        help="SnpEff annotated VCF"
    )
    parser.add_argument(
        "-o",
        "--output",
        default="variant_density.png",
        help="output figure"
    )
    parser.add_argument(
        "--window",
        type=int,
        default=10000,
        help="window size in bp (default 10000)"
    )
    parser.add_argument(
        "--detail",
        action="store_true",
        help="plot individual SnpEff effects"
    )
    args = parser.parse_args()
    positions = []
    categories = []
    effects_all = []
    with open_file(args.vcf) as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            fields=line.rstrip().split("\t")
            pos=int(fields[1])
            info=fields[7]
            effects=parse_ann(info)
            category=classify_effects(effects)
            positions.append(pos)
            categories.append(category)
            effects_all.append(effects)
    df=pd.DataFrame(
        {
            "position":positions,
            "category":categories,
            "effects":effects_all
        }
    )
    print("\nVariant summary")
    print("----------------")
    print(df.category.value_counts())
    max_position=df.position.max()
    bins=np.arange(
        0,
        max_position + args.window,
        args.window
    )
    centers=(bins[:-1]+bins[1:])/2
    plt.figure(figsize=(14,6))
    for category in [
        "gene_synonymous",
        "gene_nonsynonymous",
        "non_gene"
    ]:
        subset=df[
            df.category==category
        ]
        counts,_=np.histogram(
            subset.position,
            bins=bins
        )
        plt.plot(
            centers,
            counts,
            label=category
        )
    if args.detail:
        detail_effects=set()
        for effects in df.effects:
            detail_effects.update(effects)
        detail_effects=sorted(detail_effects)
        for effect in detail_effects:
            subset_positions=[]
            for pos,effects in zip(
                df.position,
                df.effects
            ):
                if effect in effects:
                    subset_positions.append(pos)
            if len(subset_positions)==0:
                continue
            counts,_=np.histogram(
                subset_positions,
                bins=bins
            )
            plt.plot(
                centers,
                counts,
                linestyle="--",
                alpha=0.7,
                label=effect
            )
    plt.xlabel(
        "Reference position (bp)"
    )
    plt.ylabel(
        f"Variants per {args.window:,} bp"
    )
    plt.xlim(
        0,
        max_position
    )
    plt.legend(
        fontsize=8
    )
    plt.tight_layout()
    plt.savefig(
        args.output,
        dpi=300
    )
    print(
        "\nSaved:",
        args.output
    )


if __name__=="__main__":
    main()
