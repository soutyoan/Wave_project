# Wavelets project

## Pre-requis
- python 3.0 (numpy, opencv, PIL, matplotlib.pyplot, progressbar)
- lena.jpg dans le dossier images

(Si besoin pour progressbar)

```
pip install progressbar2
```

## Structure

- doc : ensemble des documents .pdf

- src : scripts Python

	- ComplexGaborFilters.py : classe Python de la décomposition complexes en filtres de Gabor
	- LPC_computation.py : implémentation de l'indice présenté par doc/BIQA_LPC.pdf
	- detect_Blur.py : comparaison de l'indice de netteté entre une image et une filtrée par Gaussien

- images : images en niveaux de gris de test

- output : images résultats (sous-bandes, graphes...)

## Exécution des scripts
L'option --help est disponible sur tous les scripts suivants.

- LPC_computation : calcul direct de l'indice de netteté d'une image en entrée.

```console
me@machine:~$ cd src/
me@machine:~$ python3 LPC_computation -i <nom_image> [-N <nb_scales>] [-M <nb_orientations>] [-C <constante>] [-B <beta>]

Exemple :
```
python3 LPC_computation -i small_lena.jpg
```

- detect_Blur : calcul direct de l'indice de netteté d'une image en entrée et d'une image floutée par filtre Gaussien d'écart-type sigma.

```console
me@machine:~$ cd src/
me@machine:~$ python3 LPC_computation -i <nom_image> [-N <nb_scales>] [-M <nb_orientations>] [-C <constante>] [-B <beta>]

Exemple :
```
python3 detect_Blur.py -i small_lena.jpg -sig <sig_value> [-B <beta>]
```

