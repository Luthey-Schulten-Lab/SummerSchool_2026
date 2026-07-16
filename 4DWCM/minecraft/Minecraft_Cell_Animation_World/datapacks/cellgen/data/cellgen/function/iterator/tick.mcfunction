# Ran as / at the iterator

# iterator has x_current, y_current, z_current coordinates. it also has x_max, y_max, z_max coordinates


execute if score @s cg.step_cooldown matches 1.. run scoreboard players remove @s cg.step_cooldown 1

execute if score @s cg.step_cooldown matches 0 run function cellgen:iterator/step with storage iterator

