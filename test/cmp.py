import pandas as pd

for data_name in ("sim_tcr_ampli_rc", "sim_tcr_batch"):
    df1 = pd.read_csv(f"trust4_result/{data_name}/TRUST_{data_name}_report.tsv", sep="\t")
    df2 = pd.read_csv("sim_tcr.stats.tsv", sep="\t", quotechar="'")

    detected = set(df1["CDR3aa"])
    gt = set(df2["ACDR3_AA"]).union(set(df2["BCDR3_AA"]))

    tp = len(gt.intersection(detected))
    fn = len(gt.difference(detected))
    fp = len(detected.difference(gt))
    precision = tp / (tp + fp)
    accuracy = tp / (tp + fn)
    print(
        "{} -- TP: {}, FP: {}, FN: {}, Precision: {}, Accuracy: {}, F1: {}".format(
            data_name,
            tp,
            fp,
            fn,
            precision,
            accuracy,
            2 * precision * accuracy / (precision + accuracy),
        )
    )
