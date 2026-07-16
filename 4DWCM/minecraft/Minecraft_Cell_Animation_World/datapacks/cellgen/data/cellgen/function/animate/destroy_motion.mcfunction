# Can be ran as/at anywhere... or at the entity itself. Doesn't matter. The power of return run!

# The first line is just to make sure the function doesn't run if there isn't a motion to destroy. It's not needed but it just lets the player know which is nice. Actually, I removed it because this could happen at the start.
# execute unless entity @s[type=marker,tag=cg.animate,distance=0.01] run return run execute as @n[type=marker,tag=cg.animate] at @s run function cellgen:animate/destroy_motion

# If you don't feel like having the tellraw, just paste the function in.

# Max size to destroy: 110, 59, 59. Of course, if you want, 128, 64, 64 is good too.
# Assuming perfect mirroring, this means that we are destroying blocks 111 to 128 (18 inclusive) and 0 - 17 (inclusive). Start at 18, 5, 5, end at 110, 59, 59. Of course I don't really trust this so we're going larger.
fill ~0 ~0 ~0 ~31 ~31 ~31 air replace #cellgen:cell_blocks
fill ~0 ~0 ~32 ~31 ~31 ~63 air replace #cellgen:cell_blocks
fill ~0 ~32 ~0 ~31 ~63 ~31 air replace #cellgen:cell_blocks
fill ~0 ~32 ~32 ~31 ~63 ~63 air replace #cellgen:cell_blocks
fill ~32 ~0 ~0 ~63 ~31 ~31 air replace #cellgen:cell_blocks
fill ~32 ~0 ~32 ~63 ~31 ~63 air replace #cellgen:cell_blocks
fill ~32 ~32 ~0 ~63 ~63 ~31 air replace #cellgen:cell_blocks
fill ~32 ~32 ~32 ~63 ~63 ~63 air replace #cellgen:cell_blocks
fill ~64 ~0 ~0 ~95 ~31 ~31 air replace #cellgen:cell_blocks
fill ~64 ~0 ~32 ~95 ~31 ~63 air replace #cellgen:cell_blocks
fill ~64 ~32 ~0 ~95 ~63 ~31 air replace #cellgen:cell_blocks
fill ~64 ~32 ~32 ~95 ~63 ~63 air replace #cellgen:cell_blocks
fill ~96 ~0 ~0 ~127 ~31 ~31 air replace #cellgen:cell_blocks
fill ~96 ~0 ~32 ~127 ~31 ~63 air replace #cellgen:cell_blocks
fill ~96 ~32 ~0 ~127 ~63 ~31 air replace #cellgen:cell_blocks
fill ~96 ~32 ~32 ~127 ~63 ~63 air replace #cellgen:cell_blocks