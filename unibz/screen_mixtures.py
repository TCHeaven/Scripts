import csv
import statistics
from collections import defaultdict

input_file = "mixture_screen_raw.txt"
output_file = "mixture_screen_summary.tsv"

data = defaultdict(lambda: {
    "dp": [],
    "called": 0,
    "het": 0,
    "dp5": 0,
    "maf5": 0,
    "maf10": 0,
    "maf20": 0,
    "maf30": 0,
    "maf40": 0,
})

with open(input_file) as f:
    for line in f:
        sample, gt, dp, ad = line.rstrip().split("\t")

        if dp == "." or dp == "0":
            continue

        dp = int(dp)

        if gt == "./.":
            continue

        data[sample]["called"] += 1
        data[sample]["dp"].append(dp)

        if gt in ("0/1", "1/0", "0|1", "1|0"):
            data[sample]["het"] += 1

        if dp >= 5:
            data[sample]["dp5"] += 1

        if ad == ".":
            continue

        ads = ad.split(",")

        if len(ads) < 2:
            continue

        try:
            ref = int(ads[0])
            alt = sum(int(x) for x in ads[1:])
        except ValueError:
            continue

        total = ref + alt

        if total == 0:
            continue

        maf = min(ref, alt) / total

        if dp >= 5:
            if maf >= 0.05:
                data[sample]["maf5"] += 1
            if maf >= 0.10:
                data[sample]["maf10"] += 1
            if maf >= 0.20:
                data[sample]["maf20"] += 1
            if maf >= 0.30:
                data[sample]["maf30"] += 1
            if maf >= 0.40:
                data[sample]["maf40"] += 1

with open(output_file, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")

    writer.writerow([
        "Sample",
        "Called_sites",
        "Mean_DP",
        "Median_DP",
        "DP>=5",
        "Het_calls",
        "Het_fraction",
        "MAF>=5%",
        "MAF>=10%",
        "MAF>=20%",
        "MAF>=30%",
        "MAF>=40%"
    ])

    for sample in sorted(data):
        d = data[sample]

        mean_dp = statistics.mean(d["dp"]) if d["dp"] else 0
        median_dp = statistics.median(d["dp"]) if d["dp"] else 0

        het_fraction = (
            d["het"] / d["called"]
            if d["called"] > 0
            else 0
        )

        writer.writerow([
            sample,
            d["called"],
            f"{mean_dp:.2f}",
            f"{median_dp:.2f}",
            d["dp5"],
            d["het"],
            f"{het_fraction:.5f}",
            d["maf5"],
            d["maf10"],
            d["maf20"],
            d["maf30"],
            d["maf40"],
        ])

print(f"Wrote {output_file}")
