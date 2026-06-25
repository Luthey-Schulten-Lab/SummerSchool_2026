"""
Serialize merged WCM trajectories into a single ensemble (.pkl + _x/_fx .npy).

Reads the combined counts_SA_fluxes_<n>.csv files produced by merge_categorize.py
(under <main_folder>/<output_folder>/), assembles them into 3-D arrays
(species/reactions x time x replicate) with WCM_analysis.merge_trajs (parallel
across replicates), and writes the ensemble back into the first output folder.

Usage:
    python pkl.py --main_folder <DIR> --output_folders <F1> [F2 ...] [--n_cpus N]

  --main_folder     base directory (same one passed to merge_categorize.py)
  --output_folders  merged category folder(s) to serialize (default: trajs_healthy)
  --n_cpus          worker processes for the parallel read (default: 8)

The number of replicates is auto-detected from the files present.
"""

import argparse
import glob
import os
import re
import sys
from datetime import datetime

import numpy as np

import WCM_analysis


def parse_args():
    ap = argparse.ArgumentParser(description='Serialize merged WCM trajectories into a pkl ensemble.')
    ap.add_argument('-mf', '--main_folder', required=True,
                    help='base directory containing the merged output folders')
    ap.add_argument('-of', '--output_folders', nargs='+', default=['trajs_healthy'],
                    help='merged category folder(s) to serialize (default: trajs_healthy)')
    ap.add_argument('-n_cpus', '--n_cpus', type=int, default=8,
                    help='worker processes for the parallel CSV read (default: 8)')
    return ap.parse_args()


def replicate_number(path):
    return int(re.search(r'_(\d+)\.csv$', os.path.basename(path)).group(1))


def collect_traj_files(main_folder, output_folders):
    """All counts_SA_fluxes_<n>.csv across the output folders, replicate-ordered."""
    files = []
    for folder in output_folders:
        pattern = os.path.join(main_folder, folder, 'counts_SA_fluxes_*.csv')
        files += sorted(glob.glob(pattern), key=replicate_number)
    return files


def main():
    args = parse_args()
    start = datetime.now()

    traj_files = collect_traj_files(args.main_folder, args.output_folders)
    if not traj_files:
        sys.exit('No counts_SA_fluxes_*.csv found under {0}/{1}'.format(
            args.main_folder, args.output_folders))

    out_dir = os.path.join(args.main_folder, args.output_folders[0]) + os.sep
    label = 'WCMensemble_' + os.path.basename(args.main_folder.rstrip('/')) if \
        os.path.basename(args.main_folder.rstrip('/')) else 'WCMensemble'

    print('Serializing {0} replicates -> {1}{2} (started {3})'.format(
        len(traj_files), out_dir, label, start))

    # build the ensemble from the explicit file list (rep numbering = file order)
    w = WCM_analysis.WCM_ensemble()
    w.traj_files = traj_files
    w.N_reps = len(traj_files)
    w.reps = np.arange(1, w.N_reps + 1, dtype=np.int32)

    w.merge_trajs(args.n_cpus)
    w.write_merged_ensemble(out_dir, label)
    del w  # free RAM before reading back

    print('Serialized in {0}'.format(datetime.now() - start))

    # read back to confirm the ensemble loads
    w = WCM_analysis.WCM_ensemble()
    w.read_merged_ensemble(out_dir, label)
    print('Ensemble OK: {0} replicates, {1} timepoints'.format(w.N_reps, w.t))


if __name__ == '__main__':
    main()
