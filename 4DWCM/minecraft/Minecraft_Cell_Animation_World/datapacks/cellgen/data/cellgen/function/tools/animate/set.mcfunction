advancement revoke @s only cellgen:set_frame
execute if score @s cg.item_cooldown matches 1.. run return fail

scoreboard players operation @n[type=marker,tag=cg.animate] cg.motion_countdown = #animate_tick cg.rule
# execute store result score @n[type=marker,tag=cg.animate] cg.motion_tick run data get entity @s SelectedItem.components."minecraft:custom_data".cg.frame

scoreboard players operation @n[type=marker,tag=cg.animate] cg.motion_tick = @s cg.motion_tick
execute store result storage motion_tick tick int 1 run scoreboard players get @s cg.motion_tick
execute at @n[type=marker,tag=cg.animate] run function cellgen:animate/tick_macro with storage motion_tick

# Set the mainhand item's damage to 30 - the player's cg.motion_tick
# (set_damage uses a durability *fraction*; motion_tick/30 leaves damage = 30 - motion_tick)
item modify entity @s weapon.mainhand cellgen:name_motion_damage
item modify entity @s hotbar.0 cellgen:set_motion_damage

tellraw @a {"text": "Frame has been set", "color": "dark_purple"}


scoreboard players set @s cg.item_cooldown 20
