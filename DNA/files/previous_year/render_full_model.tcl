# Get the number of frames from the top molecule
set molid [molinfo top get id]
set num_frames [molinfo $molid get numframes]

# Output directory for rendered frames
set output_dir "frames"
file mkdir $output_dir

# Loop through each frame
for {set i 0} {$i < $num_frames} {incr i} {
    animate goto $i
    set frame_number [format "%04d" $i]
    set filename "$output_dir/frame$frame_number.tga"
    
    render TachyonInternal  $filename
}

