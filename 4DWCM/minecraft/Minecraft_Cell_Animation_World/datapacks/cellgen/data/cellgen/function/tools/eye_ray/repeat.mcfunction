
particle crit ~ ~ ~ 0 0 0 0 1 normal

scoreboard players add #distance scratch 1
execute if score #distance scratch matches ..200 if block ~ ~ ~ air positioned ^ ^ ^0.30 run return run function cellgen:tools/eye_ray/repeat

#playsound block.note_block.pling master @s ~ ~ ~ 1 2

return run function cellgen:tools/eye_ray/read_block_data
