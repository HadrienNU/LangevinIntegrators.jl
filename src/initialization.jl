#=
initialization.jl

To get initial state for the trajectories

Il faut créer une structure qui contient les paramètres de l'initialisation

Type d'initialisation:
    -Constant
    - Array (on fournit un tableau et ça génère successivement les CI)
    - Random:
        Uniform
        Gaussien
        Selon une pmf (si on a un potentiel et des bornes, dans ce cas ça dérive de Uniform)

En fait on a 2 types d'initialisation en série et en parallèle.
Ca veut dire qu'il y a un état interne qui détermine l'initialisation d'après
Ca sera une seed pour le cas uniforme et un index pour le tableau
Pour se simplifier la vie, on peut supposer qu'il n'y a que des initializer serial et ce sont des états qu'on passe au cas distribués.
Ca permet un meileur controle de ce qu'on passe
=#
