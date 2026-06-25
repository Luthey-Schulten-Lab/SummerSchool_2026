"""
Serialize the per-replicate 4DWCM trajectory CSVs into one memory-mappable
ensemble so the notebook loads instantly instead of re-reading + merging the
raw ~13 GB of CSVs on every run.

Reads  <data_dir>/counts_and_fluxes.<n>.csv  for n = 1..n_reps
Writes <data_dir>/<label>_WCMensemble.pkl        (small metadata: species/rxns maps,
                                                  time axis, volumes, conc factors)
       <data_dir>/<label>_WCMensemble_x.npy       (counts:  N_species x Nt x N_reps, int32)
       <data_dir>/<label>_WCMensemble_fx.npy      (fluxes:  N_rxns   x Nt x N_reps, float32)

The notebook then loads it with WCM_ensemble.read_merged_ensemble(data_dir, label),
which memory-maps the .npy files (no large RAM spike).

Usage (run from SummerSchool_2026/4DWCM/ so the analysis_scripts package imports):
    python -m analysis_scripts.serialize [data_dir] [label] [n_reps] [n_workers]
Defaults: data/  ensemble  50  8

n_workers is intentionally modest: each worker holds one ~250 MB trajectory while
loading, and all of them plus the merged array must fit in RAM.
"""

import os
import sys
import time

import numpy as np

from analysis_scripts import WCM_analysis


def main():
    data_dir = sys.argv[1] if len(sys.argv) > 1 else 'data/'
    label = sys.argv[2] if len(sys.argv) > 2 else 'ensemble'
    n_reps = int(sys.argv[3]) if len(sys.argv) > 3 else 50
    n_workers = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    if not data_dir.endswith(os.sep):
        data_dir += os.sep

    reps = np.arange(1, n_reps + 1, dtype=np.int32)
    t0 = time.time()

    w = WCM_analysis.WCM_ensemble()
    w.set_traj_files(data_dir, 'counts_and_fluxes', reps)

    print('[1/3] loading {0} trajectories ({1} workers) ...'.format(n_reps, n_workers))
    w.load_trajs(n_workers=n_workers)

    print('[2/3] merging into 3-D arrays ...')
    w.merge_trajs()

    print('[3/3] writing ensemble -> {0}{1}_WCMensemble.{{pkl,_x.npy,_fx.npy}} ...'.format(data_dir, label))
    w.write_merged_ensemble(data_dir, label)

    print('done in {0:.1f}s : N_reps={1} Nt={2} N_species={3} N_rxns={4}'.format(
        time.time() - t0, w.N_reps, w.Nt, w.N_species, w.N_rxns))


if __name__ == '__main__':
    main()
