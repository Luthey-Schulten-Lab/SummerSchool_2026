# Will be forced to run at the marker. Will generate one if one doesn't exist 
execute unless entity @s[type=marker,tag=cg.animate,distance=..0.01] run return run function cellgen:animate/check_marker

function cellgen:animate/destroy_motion

scoreboard players set @s cg.motion_tick 0
scoreboard players operation @s cg.motion_countdown = #animate_tick cg.rule
scoreboard players set @s cg.playing 1

execute store result storage motion_tick tick int 1 run scoreboard players get @n[type=marker,tag=cg.animate] cg.motion_tick

place template cellgen:animate/lm_sim0 ~ ~ ~

tag @e[type=marker,tag=cg.animate] remove init
