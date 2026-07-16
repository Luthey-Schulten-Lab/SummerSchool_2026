tellraw @a {text:"Animation: Marker already exists!", color:"#e70b0b"}

scoreboard players set @s cg.motion_tick 0
scoreboard players operation @s cg.motion_countdown = #animate_tick cg.rule
scoreboard players set @s cg.playing 1

execute store result storage motion_tick tick int 1 run scoreboard players get @n[type=marker,tag=init] cg.motion_tick

place template cellgen:animate/lm_sim0 ~ ~ ~