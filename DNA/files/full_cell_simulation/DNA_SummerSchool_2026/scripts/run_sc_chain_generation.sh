#!/bin/bash

cd "$(dirname "$0")"

scchainDirectory="/Software/"
inputDirectory="./"
outputDirectory="../data/coords/"

input_fname="${inputDirectory}Syn3A_chromosome_init.inp"
cp $input_fname $outputDirectory
log_fname="${outputDirectory}log_init.log"

SC_SEED="${SIM_SEED:-10}"
DNA_executable="${scchainDirectory}sc_chain_generation/src/gen_sc_chain --i_f=${input_fname} --o_d=${outputDirectory} --o_l=Syn3A_chromosome_init --s=${SC_SEED} --l=${log_fname} --n_t=8 --bin --xyz"

echo "Executing command: $DNA_executable"
$DNA_executable

cp ${outputDirectory}x_chain_Syn3A_chromosome_init_rep00001.bin ${outputDirectory}dna_summerschool_0.bin
cp ${outputDirectory}x_obst_Syn3A_chromosome_init_rep00001.bin ${outputDirectory}ribo_summerschool_0.bin
