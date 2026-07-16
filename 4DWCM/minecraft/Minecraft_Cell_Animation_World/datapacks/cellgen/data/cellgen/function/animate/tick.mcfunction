# Ran at / as the moving marker

execute if score @s cg.playing matches 0 run return fail


scoreboard players remove @s[scores={cg.motion_countdown=1..,cg.motion_tick=..30,cg.playing=1}] cg.motion_countdown 1

execute if score @s cg.motion_countdown matches 0 run function cellgen:animate/new_frame
