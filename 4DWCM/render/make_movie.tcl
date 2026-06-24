# =============================================================================
#  make_movie.tcl  --  4DWCM (STC-QCB summer school 2026)
#
#  Renders the currently-loaded, represented, camera-tuned scene to an image
#  sequence (one frame per synced trajectory frame) and optionally encodes it to
#  an .mp4.  Changes ONLY lighting + render quality -- it does NOT touch
#  molecules, representations, or the camera.
#
#  Target: Open OnDemand interactive desktop on a Delta A100 or A40 GPU node.
#    * Runs in the VMD GUI TkConsole (no batch -e needed).
#    * Uses the in-process OptiX ray tracer (TachyonLOptiXInternal), which runs
#      on BOTH A100 (CUDA path, no RT cores) and A40 (hardware RT cores).
#
#  Run order in the VMD TkConsole (cwd = render):
#      source load_and_sync.tcl       ;# load + sync, sets $lm_mol/$dna_mol
#      source representations.tcl      ;# visuals + camera
#      set movie_name celldiv          ;# OPTIONAL: name the output (default: mincell)
#      source make_movie.tcl          ;# THIS: lighting + render + encode
# =============================================================================

# ------------------------------- PARAMETERS ----------------------------------
# Image-sequence output directory.  Defaults to a shared, group-writable area
# (group delta_bgvl, setgid) with a PER-USER subdir so bgvl users don't clobber
# each other's frames.  Override with `setenv MOVIE_OUTDIR /path` or
# `set outdir /path` before sourcing.
if {[info exists ::env(MOVIE_OUTDIR)]} {
    set outdir $::env(MOVIE_OUTDIR)
} elseif {![info exists outdir]} {
    set _movies /projects/bgvl/SummerSchool_2026/4DWCM/render/movies
    set outdir [file join $_movies $::env(USER)]
}
# Movie name -- override by `set movie_name <name>` before sourcing (default: mincell).
# Sets the frame stem and output file, e.g. celldiv -> celldiv.00000.tga / celldiv.mp4.
if {[info exists movie_name] && $movie_name ne ""} {
    set basename $movie_name
} else {
    set basename mincell
}
set renderer   auto                    ;# auto = pick best available; or force one:
                                        ;#   TachyonLOptiXInternal  (A100 GPU)
                                        ;#   TachyonLOSPRayInternal (CPU, high quality)
                                        ;#   TachyonInternal        (CPU)
                                        ;#   snapshot               (GL grab, fastest)
set width      1920                     ;# render size; capped to the OOD desktop size
set height     1080                     ;#   for the in-process renderers (see note below)
set fps        30                       ;# playback frame rate for encoding
set stride     1                        ;# render every Nth loaded frame
set resume     1                        ;# 1 = skip frames already on disk (survive
                                        ;#     OOD disconnects/timeouts), 0 = re-render
set encode     1                        ;# 1 = encode mp4 at the end, 0 = frames only
set crf        18                        ;# x264 quality (lower = better; 18 ~ visually lossless)

# Quality look (overrides the load-time display settings for the final render).
set use_ao       1      ;# ambient occlusion (soft contact shadows) -- big quality win
set use_shadows  1      ;# hard shadows from the lights
set use_dof      0      ;# depth-of-field blur (cinematic; slower)

# Renderer candidates, best-first.  Probed at run time; the first one that
# actually writes an image is used for the whole movie.
set _candidates {TachyonLOptiXInternal TachyonLOSPRayInternal TachyonInternal snapshot}
if {$renderer ne "auto"} { set _candidates [list $renderer] }
# -----------------------------------------------------------------------------


# ------------------------------- LIGHTING ------------------------------------
# VMD has 4 directional lights (0-3).  0 = key (on by default).  Add a fill so
# the shadowed side isn't black; leave 2/3 off for a clean two-point setup.
light 0 on
light 1 on
light 2 off
light 3 off


# --------------------------- RENDER QUALITY ----------------------------------
display rendermode GLSL                 ;# AO/shadows in the GL preview + snapshot fallback
color Display Background white

if {$use_ao} {
    display ambientocclusion on
    display aoambient 0.80
    display aodirect  0.30
} else {
    display ambientocclusion off
}
display shadows [expr {$use_shadows ? "on" : "off"}]

if {$use_dof} {
    display dof on
    display dof_fnumber   64.0
    display dof_focaldist 0.70
} else {
    display dof off
}

axes location Off                       ;# no axis widget baked into the movie

# Size the window so the in-process ray tracer renders at the wanted resolution.
# NOTE: TachyonLOptiXInternal / TachyonInternal / snapshot render at the VMD GL
# window size, which is CAPPED by the OOD desktop resolution.  For true 1080p,
# make the OOD desktop session >= 1920x1080.  For resolution INDEPENDENT of the
# desktop, use the external `Tachyon` renderer with `-res W H` (offscreen) -- ask
# if you need that path.
display resize $width $height
display update


# ------------------------------- RENDER LOOP ---------------------------------
# Render to $fn and report success by checking the file was actually written
# (works whether `render` throws or silently no-ops on an unknown method).
proc _try_render {method fn} {
    catch {file delete -- $fn}
    catch {render $method $fn}
    return [expr {[file exists $fn] && [file size $fn] > 0}]
}

set last [expr {[molinfo top get numframes]-1}]

# Create the output dir and confirm WE can actually write there.  The in-process
# renderers report a misleading "Could not open file ... for writing!" (and the
# script then claims "no working renderer") when the outdir isn't writable --
# usually because VMD was launched from someone else's directory.
if {[catch {file mkdir $outdir} _e]} {
    error "make_movie: cannot create output dir '$outdir': $_e\
          \n  Fix: setenv MOVIE_OUTDIR /projects/bgvl/$::env(USER)/mincell_movie\
          \n       (or 'cd' to a directory you own), then re-source."
}
set outdir [file normalize $outdir]
set _probe [file join $outdir .write_test]
if {[catch {set _fh [open $_probe w]} _e]} {
    error "make_movie: output dir is not writable by you: $outdir\
          \n  ($_e)\
          \n  Fix: setenv MOVIE_OUTDIR /projects/bgvl/$::env(USER)/mincell_movie\
          \n       (or 'cd' to a directory you own), then re-source."
}
close $_fh
file delete -- $_probe

# Build the frame list up front (so ETA knows the total).
set frames {}
for {set f 0} {$f <= $last} {incr f $stride} { lappend frames $f }
set total [llength $frames]
puts "make_movie: $total frames (0..$last stride $stride) -> $outdir/  resume=$resume"

set active ""
set done   0
set t0     [clock seconds]
foreach f $frames {
    set fn [format "%s/%s.%05d.tga" $outdir $basename $f]

    # resume: skip frames already rendered
    if {$resume && [file exists $fn] && [file size $fn] > 0} {
        incr done
        continue
    }

    animate goto $f                     ;# step ALL molecules to frame f (they are synced)
    display update                      ;# force the scene to the current frame

    if {$active eq ""} {
        # First real render: probe candidates, keep the one that works.
        foreach m $_candidates {
            if {[_try_render $m $fn]} { set active $m; break }
        }
        if {$active eq ""} {
            error "make_movie: no working renderer among {$_candidates}"
        }
        puts "make_movie: renderer = $active"
    } elseif {![_try_render $active $fn]} {
        error "make_movie: render failed at frame $f with $active"
    }

    incr done
    if {$done == 1 || $f % 25 == 0 || $done == $total} {
        set el  [expr {[clock seconds]-$t0}]
        set eta [expr {$done>0 ? int(($total-$done)*(double($el)/$done)) : 0}]
        puts [format "  %d/%d  frame %d  elapsed %ds  eta ~%ds" $done $total $f $el $eta]
    }
}
puts "make_movie: done, $total images in $outdir/"


# --------------------------------- ENCODE ------------------------------------
# mp4 if ffmpeg is on PATH (on Delta: `module load ffmpeg` before launching VMD);
# otherwise encode with the Python imageio/Pillow modules -- part of any standard
# conda env, no ffmpeg binary needed (mp4 if imageio-ffmpeg is installed, else an
# animated GIF).  No separate shell script.
#
# glob input tolerates a stride>1 frame numbering; the scale filter forces even
# WxH, which libx264/yuv420p requires.
set _encode_py {
import sys, os, glob
outdir, base = sys.argv[1], sys.argv[2]
fps, crf, out = int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
files = sorted(glob.glob(os.path.join(outdir, base + ".*.tga")))
if not files:
    raise SystemExit("encode: no frames in %s/%s.*.tga" % (outdir, base))
ok = False
try:                                      # mp4 via imageio (uses imageio-ffmpeg if present)
    import imageio.v2 as iio
    w = iio.get_writer(out, fps=fps, macro_block_size=None)
    for f in files:
        w.append_data(iio.imread(f))
    w.close()
    print("encode: wrote %s (%d frames)" % (out, len(files)))
    ok = True
except Exception as e:
    print("encode: mp4 via imageio unavailable (%s); writing GIF instead" % type(e).__name__)
if not ok:                                # animated GIF via Pillow -- no ffmpeg needed
    if os.path.exists(out):               # drop the empty mp4 imageio left behind
        try: os.remove(out)
        except OSError: pass
    from PIL import Image
    gif, maxw, ims = os.path.splitext(out)[0] + ".gif", 960, []
    for f in files:
        im = Image.open(f).convert("RGB")
        if im.width > maxw:
            im = im.resize((maxw, int(im.height * maxw / im.width)))
        ims.append(im)
    ims[0].save(gif, save_all=True, append_images=ims[1:],
                duration=int(round(1000.0 / fps)), loop=0, optimize=True)
    print("encode: wrote %s (%d frames)" % (gif, len(files)))
}

if {$encode} {
    set out ${basename}.mp4
    set ff  [auto_execok ffmpeg]
    if {[llength $ff]} {
        puts "make_movie: encoding $out with [lindex $ff 0] ..."
        exec [lindex $ff 0] -y -framerate $fps -pattern_type glob \
             -i "$outdir/$basename.*.tga" \
             -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
             -c:v libx264 -pix_fmt yuv420p -crf $crf $out 2>@ stdout
        puts "make_movie: wrote $out"
    } else {
        set py ""
        foreach c {python python3} {
            set p [auto_execok $c]
            if {[llength $p]} { set py [lindex $p 0]; break }
        }
        if {$py eq ""} {
            puts "make_movie: no ffmpeg and no python found; frames are in $outdir/"
        } else {
            puts "make_movie: ffmpeg not on PATH; encoding via Python ($py) ..."
            exec $py -c $_encode_py $outdir $basename $fps $crf $out 2>@ stdout
        }
    }
}
