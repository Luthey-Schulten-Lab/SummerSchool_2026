advancement revoke @s only cellgen:forward
execute if score @s cg.modify_cooldown matches 1.. run return fail

execute if score @s cg.motion_tick matches ..29 run scoreboard players add @s cg.motion_tick 1

item modify entity @s hotbar.4 cellgen:name_motion_damage

scoreboard players set @s cg.modify_cooldown 5
