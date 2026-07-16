# Ran as / at the iterator

# iterator has x_current, y_current, z_current coordinates. it also has x_max, y_max, z_max coordinates

# for each value of x
#   for z
#       for y
#           perform task
#           step y
#       step z
#   step x

# Default step cooldown is 1
scoreboard players operation @s cg.step_cooldown = #default_step_cooldown cg.rule

execute store result storage iterator x int 1 run scoreboard players get @s cg.x
execute store result storage iterator y int 1 run scoreboard players get @s cg.y
execute store result storage iterator z int 1 run scoreboard players get @s cg.z

execute if entity @s[tag=cg.destroy] run place template cellgen:air
execute unless entity @s[tag=cg.destroy] run function cellgen:iterator/place_macro with storage iterator


scoreboard players add @s cg.y 1
tp @s ~ ~48 ~

# execute unless score @s cg.z >= @s cg.z_max run forceload add ~ ~48 ~48 ~96



execute if score @s cg.y > @s cg.y_max run function cellgen:iterator/large_y

# Double steps if the function has failed in the past.
# execute if score @s cg.function_passed matches 0 run function cellgen:iterator/step

#execute unless entity @p[distance=..200] run tp @p @s
