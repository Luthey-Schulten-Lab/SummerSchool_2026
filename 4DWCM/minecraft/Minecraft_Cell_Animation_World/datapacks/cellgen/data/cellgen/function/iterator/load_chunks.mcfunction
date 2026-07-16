# Note that the at has been already set
$execute positioned ~$(old_x48) ~ ~$(old_z48) run forceload remove ~ ~ ~48 ~48
$execute positioned ~$(x48) ~ ~$(z48) run forceload add ~ ~ ~48 ~48
# $tellraw @a "Forceloaded chunks ~$(x48) ~$(z48)"