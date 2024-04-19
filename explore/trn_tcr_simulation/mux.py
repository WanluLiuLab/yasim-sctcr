import pandas as pd
import os

if __name__ == "__main__":
    for t_cell_num in (100, 500, 1000):
        for replicate_num in range(1, 11):
            _sample_name = f"HU_0043_Blood_10x_sim_tcell_only_ncells{t_cell_num}_rep{replicate_num}.d"
            t_cell_stats = pd.read_csv(
                os.path.join("sim", _sample_name, "sim_t_cell.nt.fa.stats.tsv"),
                sep="\t",
            )
            (
                pd
                .read_csv(
                    os.path.join("sim", _sample_name, "sim_tcr.stats.tsv"),
                    sep="\t",
                )
                .set_index("UUID")
                .join(t_cell_stats.groupby("UUID").count(), how="inner")
                .reset_index()
                .to_csv(
                    os.path.join(
                        "sim",
                        "HU_0043_Blood_10x_sim_tcell_only_mux",
                        _sample_name+".sim_tcr.stats.withcounts.tsv"
                    ),
                    sep="\t"
                )
            )
