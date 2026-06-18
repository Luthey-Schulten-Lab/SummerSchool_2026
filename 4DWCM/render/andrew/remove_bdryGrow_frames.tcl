set currentBatch {}

# Iterate through all frames in reverse order
for {set i [molinfo top get numframes]} {$i >= 0} {incr i -1} {
    # Select anchors and hinges for the current frame
    set anchor [atomselect top "vy==7" frame $i]
    set hinge [atomselect top "vy==8" frame $i]

    # Check if the frame is part of a BdryGrow batch
    if {[$anchor num] == 0 && [$hinge num] == 0} {
        # Add frame to the current batch
        lappend currentBatch $i
    } else {
        # Process the current batch if the frame has anchors or hinges
        if {[llength $currentBatch] > 0} {
            # Retain the last frame of the batch
            set lastBdryGrow [lindex $currentBatch end]
            foreach frame $currentBatch {
                if {$frame != $lastBdryGrow} {
                    animate delete beg $frame end $frame
                }
            }
            # Clear the batch
            set currentBatch {}
        }
    }

    # Clean up selections
    $anchor delete
    $hinge delete
}

# Handle any remaining batch at the end
if {[llength $currentBatch] > 0} {
    set lastBdryGrow [lindex $currentBatch end]
    foreach frame $currentBatch {
        if {$frame != $lastBdryGrow} {
            animate delete $frame
        }
    }
}

# Update the display
display resetview

