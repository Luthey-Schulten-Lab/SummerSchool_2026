# requires x, y, z, name
kill @e[type=marker,tag=cg.iterator_marker]
kill @e[type=armor_stand,tag=cg.iterator]

scoreboard players add #id cg.id 1

$summon armor_stand ~ ~ ~ {Marker:true,NoGravity:true,Tags:["init", "cg.iterator", "cg.$(type)"],CustomName:{text:"Cell Builder", color:"#f6ab16"}}
#function cellgen:fill/breast/break with entity @s CustomName 
# This stores original point

forceload add ~ ~
summon marker ~ ~ ~ {Tags:["cg.iterator_marker","init2"]}
scoreboard players operation @n[type=marker,tag=init2] cg.id = #id cg.id
tag @n[type=marker,tag=init2] remove init2


$data modify storage iterator name set value "$(name)"

$execute as @n[type=armor_stand,tag=init] run function cellgen:iterator/set_score {"x": $(x), "y": $(y), "z": $(z)}

# Fills in the forceload
forceload add ~ ~ ~ ~

tellraw @a "Initialized Iteration"