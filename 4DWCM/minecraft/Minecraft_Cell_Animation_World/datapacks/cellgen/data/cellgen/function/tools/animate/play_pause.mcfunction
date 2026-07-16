
advancement revoke @s only cellgen:play_pause
execute if score @s cg.item_cooldown matches 1.. run return fail


# 0 = currently paused. 1 = currently playing.
execute as @n[type=marker,tag=cg.animate] at @s run function cellgen:tools/animate/marker_play_pause

scoreboard players set @s cg.item_cooldown 20
