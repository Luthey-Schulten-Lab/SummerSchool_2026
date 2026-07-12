import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.animation as animation
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Rectangle
import numpy as np
import re
import shutil
from collections import defaultdict
import ast


def get_ffmpeg_path():
    """Return ffmpeg executable path, or None if unavailable."""
    ffmpeg = shutil.which("ffmpeg")
    if ffmpeg:
        return ffmpeg

    # Common Spack install on NCSA Delta when the module is loaded.
    spack_ffmpeg = "/sw/rh9.4/spack/v1.0.0/sw/linux-x86_64_v2/ffmpeg-7.1-3avnbo4/bin/ffmpeg"
    if shutil.which(spack_ffmpeg) or __import__("os").path.isfile(spack_ffmpeg):
        return spack_ffmpeg

    return None

def parse_loops_file(filename):
    """Parses loops.txt and extracts loop start/end pairs and replication forks for animation."""
    with open(filename, 'r') as file:
        lines = file.readlines()

    frames = []
    current_loops = []
    current_forks = []

    for line in lines:
        line = line.strip()
        if line.startswith("Number of loops:"):
            if current_loops or current_forks:
                frames.append((current_loops, current_forks))
                current_loops = []
                current_forks = []
        elif line.startswith("Replication forks:"):
            try:
                fork_vals = line.split(":")[1].strip().split(",")
                current_forks = [int(f.strip()) for f in fork_vals if f.strip()]
            except ValueError:
                current_forks = []
        elif line:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    start = int(parts[0])
                    end = int(parts[1])
                    #block_start = int(parts[2]) #bug
                    #block_end = int(parts[3]) #bug
                    
                    neighbor1 = ast.literal_eval(parts[4])
                    neighbor2 = ast.literal_eval(parts[5])
                    block_start = int(neighbor1[0])      # First entry of neighbor1
                    block_end = int(neighbor2[-1])       # Last entry of neighbor2
                    current_loops.append((start, end, block_start, block_end))
                except ValueError:
                    continue

    if current_loops or current_forks:
        frames.append((current_loops, current_forks))

    every = 1
    return frames[::every]


def update(frame_data, ax, chromosome_length):
    loops, forks = frame_data
    p1 = forks[0]
    p2 = forks[1]
    p3 = 54338 + forks[0]
    p4 = 108676 - forks[0]
    ax.clear()

    ax.set_xlim(0, 108676)
    ax.set_ylim(0, 108676)
    
    if forks[0] > 24000: # get rid of forks
        ticks = [(p1+p2)//2,
        (p2+p3)//2,
        (p3+p4)//2]
        tick_labels = [r'ori$_l$'+f':\n{(p1+p2)//200}',
        r'ter:'+f'\n0',
        r'ori$_r$'+f':\n{(54338+(forks[1]-forks[0]-1)//2)//100}']
    elif forks[0] < 24000 and forks[0] > 22000: # get rid of oriL/R
        ticks = [p1,
        p2,
        (p2+p3)//2, 
        p3, 
        p4]
        tick_labels = [r'f$_{ll}$' + f':\n{p1//100}', 
        r'f$_{lr}$'+f':\n{p2//100}', 
        r'ter:'+f'\n0',
        r'f$_{rl}$'+f':\n{543}', 
        r'f$_{rr}$'+f':\n{p4//100}'] 
    elif forks[0] < 22000 and forks[0] > 6000: # have all of them
        ticks = [p1,
        (p1+p2)//2, 
        p2,
        (p2+p3)//2, 
        p3, 
        (p3+p4)//2,
        p4]
        tick_labels = [r'f$_{ll}$' + f':\n{p1//100}', 
        r'ori$_l$'+f':\n{(p1+p2)//200}',
        r'f$_{lr}$'+f':\n{p2//100}', 
        r'ter:'+f'\n0',
        r'f$_{rl}$'+f':\n{543}', 
        r'ori$_r$'+f':\n{(54338+(forks[1]-forks[0]-1)//2)//100}',
        r'f$_{rr}$'+f':\n{p4//100}']
    elif forks[0] < 6000 and forks[0] > 2000: # get rid of ter
        ticks = [p1,
        (p1+p2)//2, 
        p2,
        p3, 
        (p3+p4)//2,
        p4]
        tick_labels = [r'f$_{ll}$' + f':\n{p1//100}', 
        r'ori$_l$'+f':\n{(p1+p2)//200}',
        r'f$_{lr}$'+f':\n{p2//100}', 
        r'f$_{rl}$'+f':\n{543}', 
        r'ori$_r$'+f':\n{(54338+(forks[1]-forks[0]-1)//2)//100}',
        r'f$_{rr}$'+f':\n{p4//100}']
    else: # have outer forks and oris and ter
        ticks = [p1,
        (p1+p2)//2, 
        (p2+p3)//2, 
        (p3+p4)//2,
        p4]
        tick_labels = [r'f$_{ll}$' + f':\n{p1//100}', 
        r'ori$_l$'+f':\n{(p1+p2)//200}',
        r'ter:'+f'\n543',
        r'ori$_r$'+f':\n{(54338+(forks[1]-forks[0]-1)//2)//100}',
        r'f$_{rr}$'+f':\n{p4//100}']
    
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(tick_labels)
    ax.set_yticklabels(tick_labels)

    # Set axes (spine) line widths
    for spine in ax.spines.values():
        spine.set_linewidth(2)  # Set to desired width (e.g., 2 points)

    # Set tick line widths
    ax.tick_params(width=2)  # Tick lines
    ax.tick_params(length=6)  # Optional: longer ticks

    # Set font size of axis labels
    ax.set_xlabel("Genomic position (kb)", fontsize=14)
    ax.set_ylabel("Genomic position (kb)", fontsize=14)

    # Set font size of tick labels
    ax.tick_params(labelsize=12)
        
    # Bottom left box
    ax.add_patch(Rectangle((p1, p1), width=p2 - p1, height=p2 - p1, facecolor='green', alpha=0.07, edgecolor='none'))
    ax.add_patch(Rectangle((p1, p1), width=p2 - p1, height=p2 - p1, facecolor='none', alpha=1, edgecolor='green',lw=2))

    # Middle box
    ax.add_patch(Rectangle((p2, p2), width=p3 - p2, height=p3 - p2, facecolor='gray', alpha=0.07, edgecolor='none'))
    ax.add_patch(Rectangle((p2, p2), width=p3 - p2, height=p3 - p2, facecolor='none', alpha=1, edgecolor='gray',lw=2))

    # Top right box
    ax.add_patch(Rectangle((p3, p3), width=p4 - p3, height=p4 - p3, facecolor='magenta', alpha=0.07, edgecolor='none'))
    ax.add_patch(Rectangle((p3, p3), width=p4 - p3, height=p4 - p3, facecolor='none', alpha=1, edgecolor='magenta',lw=2))
    
    print(forks[0])
  
    # left_daughter_loops = [(s, e) for s, e in loops if s > forks[0] or e > forks[1]]
    right_daughter_loops = [(s+forks[0], e+forks[0], b1, b2) for s, e, b1, b2 in loops if s > chromosome_length or e > chromosome_length]
    mother_loops = [(forks[0]+(s-forks[0]) % chromosome_length, forks[0]+(e-forks[0]) % chromosome_length, b1, b2) for s, e, b1, b2 in loops if s <= chromosome_length and e <= chromosome_length]
    
    loops = mother_loops + right_daughter_loops
    
    red = "#d73027"
    yellow = "#e6ac00"
    green = "#1a9850"
    #red = 'red'
    #yellow = 'yellow'
    #green = 'green'
    red = '#e6194b'
    yellow =  '#ffe119'
    green ='#3cb44b'
    ms = 5.5
    
    for start, end, b1, b2 in loops:
        if b1 == 1 and b2 == 1:
            colori = red
            colorj = red
            markeri='s'
            markerj='s'
        elif b1 == 1 and b2 == 0: #start is blocked
            colori = yellow #red
            colorj = yellow
            markeri= '^'
            markerj= '>'
        elif b1 == 0 and b2 == 1: #end is blocked
            colori = yellow
            colorj = yellow #red
            markeri= '<'
            markerj='v'
        else:
            colori = green
            colorj = green
            markeri='o'
            markerj='o'
        ax.plot([start], [end], marker=markeri, color=colori, markersize=ms)
        ax.plot([end], [start], marker=markerj, color=colorj, markersize=ms)
    
def create_animation(loop_states, output_file="loops_summerschool.mp4", fps=50):
    ffmpeg_path = get_ffmpeg_path()
    if ffmpeg_path is None:
        raise RuntimeError(
            "ffmpeg is required to save .mp4 files but was not found on PATH. "
            "On Delta, run: module load ffmpeg/7.1"
        )

    mpl.rcParams["animation.ffmpeg_path"] = ffmpeg_path
    if not animation.writers.is_available("ffmpeg"):
        raise RuntimeError(f"matplotlib cannot use ffmpeg at {ffmpeg_path}")

    fig, ax = plt.subplots(figsize=(10, 10))
    ani = animation.FuncAnimation(
        fig,
        update,
        frames=loop_states,
        fargs=(ax, 54338),
        repeat=True
    )
    writer = animation.FFMpegWriter(fps=fps, codec="libx264", extra_args=["-pix_fmt", "yuv420p"])
    print(f"Saving {len(loop_states)} frames to {output_file} ...")
    ani.save(output_file, writer=writer)
    print(f"Saved {output_file}")
    plt.show()


if __name__ == "__main__":
    loops_data = parse_loops_file("./DNA_SummerSchool_2026/data/loops/loops_summerschool.txt")
    create_animation(loops_data)
