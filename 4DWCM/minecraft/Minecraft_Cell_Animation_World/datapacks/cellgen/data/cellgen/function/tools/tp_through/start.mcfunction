#Ran As, at, rotated as
advancement revoke @s only cellgen:tp_through

execute if score @s cg.item_cooldown matches 1.. run return fail
scoreboard players set @s cg.item_cooldown 15


scoreboard players set #distance scratch 0
execute anchored eyes rotated ~ ~ positioned ^ ^ ^ run function cellgen:tools/tp_through/repeat
execute at @s run playsound entity.enderman.teleport master @s ~ ~ ~ 1 1.9