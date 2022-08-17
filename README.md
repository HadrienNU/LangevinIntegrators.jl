# LangevinIntegrators

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://HadrienNU.github.io/LangevinIntegrators.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://HadrienNU.github.io/LangevinIntegrators.jl/dev/)
[![Build Status](https://github.com/HadrienNU/LangevinIntegrators.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/HadrienNU/LangevinIntegrators.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/HadrienNU/LangevinIntegrators.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/HadrienNU/LangevinIntegrators.jl)

Dossier principal pour les sources de l'intégrateur stochastique en julia



Code propre d'intégration des équations de Langevin en C++ ou Julia. (Hadrien)
    
    S'il reste du temps, on pourra regarder l'optimisation d'un code pour la génération de trajectoire à partir d'une équation de Langevin.
    
    Besoin :
    Le code doit générer des trajectoires à partir d'une équation de Langevin dont les paramètres sont fournis par un fichier extérieur. Idéalement le code doit être modulable pour pouvoir être modifié plus tard pour inclure plus de features.
    
    
    Implémentation de diverses bases functionnelles, mais idéalement aussi une équation mathématiques arbitraire.
    Au choix, overdamped ou underdamped Langevin.
    Un code qui parallélise la génération de trajectoire (de manière simple, une trajectoire par processus) (ou bien une seule trajectoire et parallèlisation avec parallel?).
    Pour avoir un code modulable, f+= new force, d*=new diff
    
    Structure du fichier de paramètres:
    
    Configurational space (boite, tore, sphere)
    Boite (nombre de dimension, conditions aux bords, overdamped/underdamped)
    
    Initial configuration: fixé ou aléatoire.
    
    Sampling : Nombre de trajectoires, nombre de pas de temps, pas de temps, choix de l'intégrateur, choix du format de sortie.
    Bonus: critère d'arrêt pour du sampling de FPT.
    
    
    Force (soit nulle, une expression functionnelle (idéalement depuis un potentiel), ou une base functionnelle avec parmètres pris d'un fichier. Quelques cas tests directement implémentés( double well potential essentiellement)? 
    
    Bain :
    Diffusion et temperature.
    Bonus: Possibilité d'avoir des bains non-markoviens dans cette section (soit par variables cachées, soit par bruit corrélés)
    

    Fourni: Code Python existant pour comparaison. https://ndsimulator.rtfd.io/
