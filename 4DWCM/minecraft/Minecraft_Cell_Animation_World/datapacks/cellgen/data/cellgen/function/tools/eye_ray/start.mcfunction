#Ran As, at, rotated as
advancement revoke @s only cellgen:info_get

execute if score @s cg.item_cooldown matches 1.. run return fail
scoreboard players set @s cg.item_cooldown 15


scoreboard players set #distance scratch 0
execute anchored eyes rotated ~ ~ positioned ^ ^ ^ run function cellgen:tools/eye_ray/repeat
