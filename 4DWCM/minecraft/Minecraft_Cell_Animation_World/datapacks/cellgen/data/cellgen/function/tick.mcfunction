execute as @e[type=armor_stand,tag=cg.iterator] at @s run function cellgen:iterator/tick
execute as @e[type=marker,tag=cg.animate] at @s run function cellgen:animate/tick

scoreboard players remove @a[scores={cg.item_cooldown=1..}] cg.item_cooldown 1
scoreboard players remove @a[scores={cg.modify_cooldown=1..}] cg.modify_cooldown 1
