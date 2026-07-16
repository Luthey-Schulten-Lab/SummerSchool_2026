
particle dust{color:16777215,scale:0.8} ~ ~ ~ 0 0 0 0 1 normal

scoreboard players add #distance scratch 1
execute if score #distance scratch matches ..50 if block ~ ~ ~ air positioned ^ ^ ^0.5 run return run function cellgen:tools/tp_through/repeat

# If above ever ends, (block not air)
scoreboard players set #distance scratch 0
function cellgen:tools/tp_through/repeat2