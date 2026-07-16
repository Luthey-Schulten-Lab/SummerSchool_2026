# We are given that this is when the y is large enough to cause problems.

# Store current location 
execute store result storage iterator old_x48 int 48 run scoreboard players get @s cg.x
execute store result storage iterator old_z48 int 48 run scoreboard players get @s cg.z
execute if score @s cg.x matches 0 if score @s cg.z matches 0 run forceload add ~ ~

# Start modifying scores
scoreboard players add @s cg.z 1
tp @s ~ ~ ~48
data modify entity @s Pos[1] set from entity @n[type=marker,tag=cg.iterator_marker] Pos[1]
scoreboard players set @s cg.y 0

# If z too large...
execute if score @s cg.z > @s cg.z_max run scoreboard players add @s cg.x 1

# If x too large don't even think about it: # Kill if x too large.
execute if score @s cg.x > @s cg.x_max run say "Finished Building!"
execute if score @s cg.x > @s cg.x_max run data modify entity @s Pos set from entity @n[type=marker,tag=cg.iterator_marker] Pos
execute if score @s cg.x > @s cg.x_max run kill @e[type=marker,tag=cg.iterator_marker]
execute if score @s cg.x > @s cg.x_max run return run kill @s




execute if score @s cg.z > @s cg.z_max run tp @s ~48 ~ ~
execute if score @s cg.z > @s cg.z_max run data modify entity @s Pos[1] set from entity @n[type=marker,tag=cg.iterator_marker] Pos[1]
# Because modifying position doesn't change where the command was run:
execute if score @s cg.z > @s cg.z_max run data modify entity @s Pos[2] set from entity @n[type=marker,tag=cg.iterator_marker] Pos[2]
execute if score @s cg.z > @s cg.z_max run scoreboard players set @s cg.z 0

# Now that the X and Z have been set, Now would be a good time to add a function that can store the data location.
execute store result storage iterator x48 int 48 run scoreboard players get @s cg.x
execute store result storage iterator z48 int 48 run scoreboard players get @s cg.z
execute at @n[type=marker,tag=cg.iterator_marker] run function cellgen:iterator/load_chunks with storage iterator


