# Programmation Parallèle

Rapport d'Arnaud Grégoire

## TD/TP1 Ensemble de Mandelbrot

### Ensemble de MandelBrot

En mathématiques, "l'ensemble de Mandelbrot est une fractale définie comme l'ensemble des points c du plan complexe pour lesquels la suite de nombres complexes définie par récurrence par :

{ z 0 = 0 z n + 1 = z n 2 + c {\displaystyle {\begin{cases}z_{0}=0\\z_{n+1}=z_{n}^{2}+c\end{cases}}} \begin{cases} z_0=0\\ z_{n+1}=z_n^2+c \end{cases}

est bornée.

L'ensemble de Mandelbrot a été découvert par Gaston Julia et Pierre Fatou1 avant la Première Guerre mondiale. Sa définition et son nom actuel sont dus à Adrien Douady, en hommage aux représentations qu'en a réalisées Benoît Mandelbrot dans les années 1980. Cet ensemble permet d'indicer les ensembles de Julia : à chaque point du plan complexe correspond un ensemble de Julia différent. Les points de l'ensemble de Mandelbrot correspondent précisément aux ensembles de Julia connexes, et ceux en dehors correspondent aux ensembles de Julia non connexes. Cet ensemble est donc intimement lié à l'ensemble de Julia, ils produisent d'ailleurs des formes similairement complexes.

Les images de l'ensemble de Mandelbrot sont réalisées en parcourant les nombres complexes sur une région carrée du plan complexe et en déterminant pour chacun d'eux si le résultat tend vers l'infini ou pas lorsqu'on y itère une opération mathématique. On considère la partie réelle et imaginaire de chaque nombre complexe comme des coordonnées et chaque pixel est coloré selon la rapidité de divergence, ou si elle ne diverge pas.

Les images de l'ensemble de Mandelbrot exposent une limite élaborée qui révèle progressivement des détails récursifs toujours plus fins en augmentant le grossissement. La limite de l'ensemble est constituée de plus petites versions de la forme principale, donc la propriété fractale de l'autosimilarité s'applique à l'ensemble tout entier (et pas simplement à certaines parties)." Wikipedia.

### Compilation

Pour compiler le code, on fait la commande

```sh
make mandel
```

Cette commande appelle un makefile qui fera la commande suivante 

```sh
mpicc mandel.c -o mandel -lm
```

### Exécution

Pour éxécuter le code, on réalise la commande :

```sh
mpirun -np [nb_processus] mandel [paramètres]
```

### Calcul séquentiel


#### Choix techniques

On représente l'image comme un tableau à 1 dimension de longeur w\*h (avec w et h la largeur et la hauteur en pixels de l'image). Etant donné que l'on stocke chaque valeur de pixel sur 1 octet (type char), on alloue donc une mémoire de w\*h\*sizeof(char).

Pour parcourir ce tableau, et afin de rendre compte de l'aspect spatial de l'ensemble de mandelbrot, on réalise 2 boucles imbriquées sur la largeur et la hauteur. Le calcul est réalisé à chaque pixel de l'image en fonction des coordonnées x et y de ce pixel, calculées à l'aide d'un incément à partir des coordonnées minimales. L'incrément dépend des limites de domaine (x mi et max, y min et max) et de la résolution de l'image.

Enfin, le paramètre profondeur correspond à la profondeur de calcul de la valeur de l'ensemble de mandelbrot à chaque pixel.



### Calcul parallèle

#### Difficulté à parallèliser

On a dit plus haut qu'il fallait optimiser les étapes dépendant de la résolution et de la profondeur. Pour cette dernière il s'agit de la fonction qui calcule la valeur de l'ensemble de mandelbrot à une position donnée. Cependant, chaque itération de cette fonction dépend de la précédente. On ne peut donc pas la parralèliser.

On se contentera donc de parallèliser le parcours du tableau correspondant à l'image, dont les itérations sont totalement indépendantes.

### Algorithme 1 : Calcul de l'ensemble du Mandelbrot en charge statique-Maitre

Paramètres : 

 - my_rank : rang du processeur
 - H_loc : longeuer de l'image locale
 - H : longueur de l'image
 - Xmin
 - Xinc
 - Ymin_loc
 - Yinc

Sortie : Calcul traité par le maître avec positionnement du pointeur au début de son image locale et allocation dynamique de l'image globale.

```c

```

## Résultats Obtenus

Commande lancé : 

```sh

mpirun -np 2 mandel 800 800 -0.736 -0.184 -0.735 -0.183 500

0.735 -0.183 500
Domaine: {[-0.736,-0.184]x[-0.735,-0.183]}
Increment : 1.25156e-06 1.25156e-06
Prof: 500
Dim image: 800x800
Domaine: {[-0.736,-0.184]x[-0.735,-0.183]}
Increment : 1.25156e-06 1.25156e-06
Prof: 500
Dim image: 800x800
1.20747
Temps total de calcul : 1.20747 sec
Temps total de calcul : 1.21253 sec
1.21253
```

![mandel](img/mandel0.ras)

Si j'augmente les dimensions de l'image 800x800 => 2500x2500 on a :

```sh
mpirun -np 2 mandel 2500 2500 -0.736 -0.184 -0.735 -0.183 500

 -0.735 -0.183 500
Domaine: {[-0.736,-0.184]x[-0.735,-0.183]}
Increment : 4.0016e-07 4.0016e-07
Prof: 500
Dim image: 2500x2500
Domaine: {[-0.736,-0.184]x[-0.735,-0.183]}
Increment : 4.0016e-07 4.0016e-07
Prof: 500
Dim image: 2500x2500
11.7508
Temps total de calcul : 11.7508 sec
Temps total de calcul : 11.7655 sec
11.7655
```
![mandel high res](img/mandel1.ras)
