#!/usr/bin/env bash
set -ue

cd sim/HU_0043_Blood_10x.sim.d
mkdir -p mux_diff_depth mux_diff_rlen
echo "${PWD}/mux_diff_depth" "${PWD}/mux_diff_rlen"

for rlen in 50 100 150 250; do
    mkdir -p mux_diff_rlen/"${rlen}"
    while read -r line; do
        cat \
            art_gex_diff_rlen/sim_gex_rlen"${rlen}"/"${line}".d_art.fq \
            > mux_diff_rlen/"${rlen}"/"${line}".fq
    done < n_cell_bc.txt 
    while read -r line; do
        cat \
            art_gex_diff_rlen/sim_gex_rlen"${rlen}"/"${line}".d_art.fq \
            art_diff_rlen/sim_tcr_rlen"${rlen}"_tcrd400.d/"${line}"_A.fq \
            art_diff_rlen/sim_tcr_rlen"${rlen}"_tcrd400.d/"${line}"_B.fq \
            > mux_diff_rlen/"${rlen}"/"${line}".fq
    done < t_cell_bc.txt
done

for depth in 2 4 8 10 20 40 60 80 100; do
    mkdir -p mux_diff_depth/"${depth}"
    while read -r line; do
        cat \
            art_gex_diff_rlen/sim_gex_rlen150/"${line}".d_art.fq \
            > mux_diff_depth/"${depth}"/"${line}".fq
    done < n_cell_bc.txt
    while read -r line; do
        cat \
            art_gex_diff_rlen/sim_gex_rlen150/"${line}".d_art.fq \
            art_diff_depth/sim_tcr_rlen150_tcrd"${depth}".d/"${line}"_A.fq \
            art_diff_depth/sim_tcr_rlen150_tcrd"${depth}".d/"${line}"_B.fq \
            > mux_diff_depth/"${depth}"/"${line}".fq
    done < t_cell_bc.txt
done

cd ..
for dir in HU_0043_Blood_10x_sim_tcell_only_*.d; do
    cd "${dir}"
    mkdir -p mux
    echo "${PWD}/mux"
    while read -r line; do
        cat \
            art_sim_gex_rlen150.d/"${line}".d_art.fq \
            art_sim_t_cell_rlen150.d/"${line}"_A.fq \
            art_sim_t_cell_rlen150.d/"${line}"_B.fq \
            > mux/"${line}".fq
    done < t_cell_bc.txt
    cd ..
done
