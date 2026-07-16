particle dust{color:0,scale:0.9} ~ ~ ~ 0 0 0 0 1 normal

scoreboard players add #distance scratch 1
execute if score #distance scratch matches ..25 unless block ~ ~ ~ air positioned ^ ^ ^1 run return run function cellgen:tools/tp_through/repeat2

tp @s ^ ^ ^0.5