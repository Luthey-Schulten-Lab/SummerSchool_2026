
scoreboard players set @s cg.x 0
scoreboard players set @s cg.y 0
scoreboard players set @s cg.z 0

scoreboard players set @s cg.step_cooldown 3

$scoreboard players set @s cg.x_max $(x)
$scoreboard players set @s cg.y_max $(y)
$scoreboard players set @s cg.z_max $(z)

scoreboard players operation @s cg.id = #id cg.id

effect give @s glowing infinite

tag @s remove init

tellraw @a "Ran successfully"