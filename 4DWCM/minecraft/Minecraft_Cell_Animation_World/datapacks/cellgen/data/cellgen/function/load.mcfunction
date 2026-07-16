tellraw @a {"text":"Cell Generator Data Pack Loaded", "bold":true, "color": "#18e230"} 
gamerule respawn_radius 0


scoreboard objectives add cg.rule dummy

# # Ensure spawn chunks are loaded
# forceload add 0 0

scoreboard objectives add cg.x dummy
scoreboard objectives add cg.y dummy
scoreboard objectives add cg.z dummy
scoreboard objectives add cg.x_max dummy
scoreboard objectives add cg.y_max dummy
scoreboard objectives add cg.z_max dummy
scoreboard objectives add cg.id dummy
scoreboard objectives add cg.step_cooldown dummy
scoreboard players set #default_step_cooldown cg.rule 2
scoreboard players set #animate_tick cg.rule 20


scoreboard objectives add cg.function_passed dummy

# A player's motion_tick counts where the player's star is
scoreboard objectives add cg.motion_tick dummy
scoreboard objectives add cg.motion_countdown dummy

scoreboard objectives add cg.playing dummy

scoreboard objectives add cg.item_cooldown dummy
scoreboard objectives add cg.modify_cooldown dummy

scoreboard objectives add scratch dummy


# Check to see if this function was ran once upon world generaton:
execute if score $check_structure cg.load matches 1 run return run tellraw @a {"text":"World Already Loaded", "bold":true, "color": "#fa5b1c"} 


scoreboard objectives add cg.load dummy
tellraw @a {"text":"Loading Structures:", "bold":true, "color": "#656565"} 

# Place your structure at spawn
setworldspawn 0 1 0
execute positioned 0 0 0 run place template cellgen:setup ~-2 ~-1 ~-5

# This code sets $check_structure = 1 no matter what might happen
scoreboard players set $check_structure cg.load 1

# WARNING: This code might not even place at all if there is nothing labeled 0 0 0.
# execute store success score $check_structure cg.load run place template cellgen:cell ~2 ~1 ~2

# To be honest, the xyz values should be a placeholder, 
# but I set it permanent because I do want this structure to be bounded.
# Picking to leave this line at the end so I can add it later.
# execute store success score $check_structure cg.load run function cellgen:iterator/start {"name": "cell/martini3_", "x": 8, "y": 8, "z": 8, "type": "build"}
