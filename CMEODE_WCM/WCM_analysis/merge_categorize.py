"""
Merge + categorize Whole-Cell Model trajectories.

For every replicate found in the given output folders (each folder is the result
of one mpirun launch and holds counts_<n>.csv / SA_<n>.csv / Flux_<n>.csv /
log_<n>.txt), this script:

  1. classifies the replicate as healthy / other / short (see METRIC below),
  2. concatenates its counts + SA + flux into one wide counts_SA_fluxes_<n>.csv,
  3. writes it (and a copy of its log) into <main_folder>/trajs_<category>/.

Replicates are renumbered 1..N per category across all input folders.

Usage:
    python merge_categorize.py --main_folder <DIR> --output_folders <F1> [F2 ...]

  --main_folder     base directory that contains the output folders
  --output_folders  one or more result folders (relative to main_folder),
                    one per mpirun launch

Reps/timepoints are auto-detected; no other arguments are needed.
"""

import argparse
import glob
import os
import re
import shutil
import sys
from datetime import datetime

import numpy as np
import pandas as pd

# Classification metric (kept as a constant so the CLI stays at two arguments):
#   0: full length, ATP & GTP accumulation == 0, pep count > 1
#   1: full length, ATP & GTP shortage    == 0, pep count > 1   (production default)
#   2: ATP accumulation >= 0 up to the surface-area doubling time
#   3: pep count > 0 up to the surface-area doubling time
# Metrics 0 and 1 also split off "short" (unfinished) replicates; 2 and 3 do not.
METRIC = 1

CATEGORIES = ('healthy', 'other', 'short')
TRAJ_DIR_PREFIX = 'trajs_'  # category dirs are <main_folder>/trajs_<category>/


def parse_args():
    ap = argparse.ArgumentParser(description='Merge + categorize WCM trajectories.')
    ap.add_argument('-mf', '--main_folder', required=True,
                    help='base directory containing the output folders')
    ap.add_argument('-of', '--output_folders', required=True, nargs='+',
                    help='result folders (one per mpirun launch), relative to main_folder')
    return ap.parse_args()


def replicate_number(path):
    """Trailing integer of e.g. counts_12.csv -> 12."""
    return int(re.search(r'_(\d+)\.csv$', os.path.basename(path)).group(1))


def header_width(csv_path):
    """Number of columns in a CSV (cheap: reads only the first line)."""
    with open(csv_path) as f:
        return len(f.readline().split(','))


def add_initial_flux(flux_path):
    """No flux is recorded at t=0, so duplicate the t=1 column into a t=0 column."""
    flux = pd.read_csv(flux_path)
    flux.insert(1, '0.0', flux[flux.columns[1]])
    return flux


def surface_area_doubling_index(surface_area, t):
    """First time index where the surface area has doubled (-1 if it never does)."""
    scaled = surface_area / surface_area[0]
    hit = np.where(scaled > 2)[0]
    return hit[0] if hit.size > 0 else -1


def _row(df, name):
    """The trajectory of a single species/row (drops the leading name column)."""
    return df.loc[df['Time'] == name].to_numpy()[0][1:]


def classify(counts, sa, full_width):
    """Return 'healthy' | 'other' | 'short' for one replicate, per METRIC."""
    finished = counts.shape[1] == full_width
    pep = _row(counts, 'M_pep_c')

    if METRIC == 0:
        atp, gtp = _row(counts, 'M_atp_c_accumulative'), _row(counts, 'M_gtp_c_accumulative')
    else:
        atp, gtp = _row(counts, 'M_atp_c_shortage'), _row(counts, 'M_gtp_c_shortage')

    if METRIC in (0, 1):
        if not finished:
            return 'short'
        if METRIC == 0:
            healthy = np.min(atp) == 0 and np.min(gtp) == 0 and np.min(pep) > 1
        else:
            healthy = np.max(atp) == 0 and np.max(gtp) == 0 and np.min(pep) > 1
        return 'healthy' if healthy else 'other'

    # metrics 2, 3: judge only up to the surface-area doubling time
    t = [int(float(c)) for c in counts.columns[1:]]
    idx = surface_area_doubling_index(_row(sa, 'SA_nm2'), t)
    if METRIC == 2:
        return 'healthy' if np.min(atp[:idx]) >= 0 else 'other'
    return 'healthy' if np.min(pep[:idx]) > 0 else 'other'


def collect_replicates(main_folder, output_folders):
    """List every replicate (counts/SA/flux/log paths) across the output folders."""
    reps = []
    for folder in output_folders:
        folder_dir = os.path.join(main_folder, folder)
        counts_files = sorted(glob.glob(os.path.join(folder_dir, 'counts_*.csv')),
                              key=replicate_number)
        if not counts_files:
            print('WARNING: no counts_*.csv in {0}'.format(folder_dir))
        for counts_path in counts_files:
            n = replicate_number(counts_path)
            reps.append({
                'folder': folder,
                'n': n,
                'counts': counts_path,
                'sa': os.path.join(folder_dir, 'SA_{0}.csv'.format(n)),
                'flux': os.path.join(folder_dir, 'Flux_{0}.csv'.format(n)),
                'log': os.path.join(folder_dir, 'log_{0}.txt'.format(n)),
            })
    return reps


def _max_existing_index(cat_dir):
    """Highest counts_SA_fluxes_<n>.csv index already in a category dir (0 if none).

    Lets a later run on a NEW output folder append (continue the numbering) instead
    of restarting at 1 and overwriting replicates already merged into that dir.
    """
    nums = [replicate_number(f)
            for f in glob.glob(os.path.join(cat_dir, 'counts_SA_fluxes_*.csv'))]
    return max(nums) if nums else 0


def main():
    args = parse_args()
    main_folder = args.main_folder
    start = datetime.now()
    print('Merging trajectories from {0} (metric {1}) at {2}'.format(
        args.output_folders, METRIC, start))

    reps = collect_replicates(main_folder, args.output_folders)
    if not reps:
        sys.exit('No replicates found under {0}'.format(main_folder))

    # a replicate is "finished" if its counts file has the full number of columns
    full_width = max(header_width(r['counts']) for r in reps)

    out_dirs = {cat: os.path.join(main_folder, TRAJ_DIR_PREFIX + cat) for cat in CATEGORIES}
    for path in out_dirs.values():
        os.makedirs(path, exist_ok=True)

    # continue numbering from replicates already merged into each category dir, so
    # running on a new output folder appends instead of overwriting (restarting at 1)
    start_counts = {cat: _max_existing_index(out_dirs[cat]) for cat in CATEGORIES}
    counters = dict(start_counts)
    for r in reps:
        counts = pd.read_csv(r['counts'])
        sa = pd.read_csv(r['sa'])
        flux = add_initial_flux(r['flux'])

        category = classify(counts, sa, full_width)
        counters[category] += 1
        idx = counters[category]

        merged = pd.concat([counts, sa, flux], ignore_index=True)
        merged.to_csv(os.path.join(out_dirs[category],
                                   'counts_SA_fluxes_{0}.csv'.format(idx)), index=False)
        if os.path.exists(r['log']):
            shutil.copyfile(r['log'], os.path.join(out_dirs[category],
                                                   'log_{0}.txt'.format(idx)))

        print('  {0}/rep {1:<3d} -> {2:<7s} #{3:<3d}  (shape {4})'.format(
            r['folder'], r['n'], category, idx, merged.shape))

    added = {cat: counters[cat] - start_counts[cat] for cat in CATEGORIES}
    print('\nClassified {0} replicates this run: {1}'.format(len(reps), added))
    if any(start_counts.values()):
        print('(appended to existing; cumulative totals per category: {0})'.format(counters))
    if len(reps):
        print('healthy fraction this run: {0:.2f}'.format(added['healthy'] / len(reps)))
    print('Done in {0}'.format(datetime.now() - start))


if __name__ == '__main__':
    main()
