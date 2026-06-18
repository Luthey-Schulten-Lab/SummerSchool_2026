# Get the total number of frames in the top molecule
set numFrames [molinfo top get numframes]

# Loop through all frames in reverse order
for {set i [expr {$numFrames - 1}]} {$i >= 0} {incr i -1} {
    # Check if the frame index is not the first or a multiple of 4
    if {$i != 0 && $i % 4 != 0} {
        animate delete beg $i end $i
    }
}

# Update the display
display resetview

