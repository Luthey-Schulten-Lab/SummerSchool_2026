#loot give @p mine ~ ~ ~ netherite_pickaxe[enchantments={silk_touch:1}]

# Originally note_block.pling with volume 1 pitch 2
execute unless block ~ ~ ~ air at @s run playsound block.note_block.hat master @s ~ ~ ~ 1 1

execute if block ~ ~ ~ glass run dialog show @s cellgen:cell/outer_cytoplasm
execute if block ~ ~ ~ red_stained_glass run dialog show @s cellgen:cell/ribosome
execute if block ~ ~ ~ red_concrete run dialog show @s cellgen:cell/ribosome
execute if block ~ ~ ~ yellow_concrete run dialog show @s cellgen:cell/dna
execute if block ~ ~ ~ lime_stained_glass run dialog show @s cellgen:cell/membrane

tellraw @a "Executed read_block_data"