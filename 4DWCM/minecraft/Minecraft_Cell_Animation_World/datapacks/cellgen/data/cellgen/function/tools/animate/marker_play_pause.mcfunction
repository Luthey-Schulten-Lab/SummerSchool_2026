execute if score @s cg.playing matches 0 if score @s cg.motion_tick matches 30.. run scoreboard players set @s cg.motion_tick -1

execute if score @s cg.playing matches 0 run scoreboard players set @s cg.playing 100
execute if score @s cg.playing matches 1 run scoreboard players set @s cg.playing 101

execute if score @s cg.playing matches 100 run scoreboard players operation @s cg.motion_countdown = #animate_tick cg.rule

execute if score @s cg.playing matches 100 run scoreboard players set @s cg.playing 1
execute if score @s cg.playing matches 101 run scoreboard players set @s cg.playing 0

execute if score @s cg.playing matches 0 run tellraw @a {"text": "Animation Paused", "color":"#f59221"}
execute if score @s cg.playing matches 1 run tellraw @a {"text": "Animation Playing", "color":"#f59221"}