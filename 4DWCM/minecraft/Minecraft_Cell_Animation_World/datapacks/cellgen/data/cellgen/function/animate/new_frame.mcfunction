scoreboard players operation @s cg.motion_countdown = #animate_tick cg.rule

# advance tick and store in storage.motion_tick.tick
scoreboard players add @s cg.motion_tick 1
execute store result storage motion_tick tick int 1 run scoreboard players get @s cg.motion_tick

# Will be run once at 30.
execute if score @s cg.motion_tick matches 31.. run return run scoreboard players set @s cg.playing 0

execute at @a[nbt={Inventory:[{Slot:0b,components:{"minecraft:custom_data":{cg.play_pause:1b}}}]}] run item modify entity @p hotbar.0 cellgen:set_motion_damage

# call the correct frame function (20 ticks per frame)



execute at @s run function cellgen:animate/tick_macro with storage motion_tick


# execute if score @s cg.motion_tick matches 30.. run funcstion cellgen:animate/destroy_motion
