advancement revoke @s only cellgen:set_speed
execute if score @s cg.item_cooldown matches 1.. run return fail

scoreboard players add #animate_tick cg.rule 10
execute if score #animate_tick cg.rule matches 70.. run scoreboard players set #animate_tick cg.rule 10
item modify entity @s weapon.mainhand cellgen:edit_speed
tellraw @a [{"text": "Ticks per Frame has been set to ", color:"#f5a742"}, { "score": { "name": "#animate_tick", "objective": "cg.rule" }, "italic": false },]

scoreboard players set @s cg.item_cooldown 15