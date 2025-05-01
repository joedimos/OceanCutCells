# Ocean CC Project

## Description
Le projet Ocean CC est une implémentation de modélisation océanique en Julia. Il comprend des fonctionnalités pour le calcul des flux, la gestion des cellules coupées et des fonctions de forçage, permettant une simulation précise des dynamiques océaniques.

## Structure du projet
- `src/ocean_cc.jl`: Code principal du projet.
- `test/runtests.jl`: Fichier de test pour vérifier le bon fonctionnement du code.
- `Project.toml`: Fichier de configuration du projet, spécifiant les dépendances.
- `Manifest.toml`: Liste complète des dépendances du projet, générée automatiquement.
- `README.md`: Documentation du projet.

## Installation
Pour installer le projet, clonez le dépôt et utilisez le gestionnaire de paquets Julia pour ajouter les dépendances spécifiées dans `Project.toml`.

```julia
using Pkg
Pkg.activate("path/to/ocean_cc_project")
Pkg.instantiate()
```

## Exécution
Pour exécuter le projet, utilisez le fichier `src/ocean_cc.jl`. Assurez-vous que toutes les dépendances sont installées et que l'environnement est activé.

## Tests
Les tests peuvent être exécutés en utilisant le fichier `test/runtests.jl`. Cela garantira que toutes les fonctionnalités du projet fonctionnent comme prévu.

```julia
include("test/runtests.jl")
```

## Contribuer
Les contributions sont les bienvenues ! Veuillez soumettre une demande de tirage avec des modifications ou des améliorations.

## License
Ce projet est sous licence MIT. Veuillez consulter le fichier LICENSE pour plus de détails.# Ocean_cc
Oceananigans 'test' integration of cut-cells
