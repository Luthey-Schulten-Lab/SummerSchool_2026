# tellraw @a {text:"Had to check marker"}

execute if entity @e[type=marker,tag=cg.animate] as @n[type=marker,tag=cg.animate] at @s run return run function cellgen:animate/marker_exists

tellraw @a [{text:"No Marker Exists: ", "color":red}, {"text": "New Marker Generated", "color": "#33c706"}]

kill @e[type=marker,tag=cg.animate]

summon marker ~ ~ ~ {Tags:[cg.animate]}
forceload add ~ ~

execute as @n[type=marker,tag=cg.animate] at @s run function cellgen:animate/start