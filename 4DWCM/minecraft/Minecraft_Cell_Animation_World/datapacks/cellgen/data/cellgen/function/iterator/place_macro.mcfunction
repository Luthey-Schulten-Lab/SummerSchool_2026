# $execute store success score @s cg.function_passed run place template cellgen:$(name)$(x)_$(y)_$(z)

execute if entity @s[tag=cg.destroy] run return run place template cellgen:air
$execute if entity @s[tag=cg.build] run place template cellgen:$(name)$(x)_$(y)_$(z)


