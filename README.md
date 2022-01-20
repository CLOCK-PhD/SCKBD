# Single Cell K-mer Bulk Decomposition

## But du projet

Le *scRNAseq* est une technique puissante qui permet de **quantifier l'abondance des transcrits au niveau cellulaire**, utile pour identifier les sous-populations cellulaires dans un échantillon. Étant plus complexe et cher que le *RNAseq* classique (*bulk RNAseq*), la quantité de données scRNAseq est beaucoup moins importante que pour le bulk.

De fait, les techniques basées sur l'expression des gènes ([Bisque](https://www.nature.com/articles/s41467-020-15816-6), [SCDC](https://pubmed.ncbi.nlm.nih.gov/31925417/), [Music](https://www.nature.com/articles/s41467-018-08023-x), [CIBERSORT](https://www.cell.com/cell-systems/fulltext/S2405-4712(16)30266-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2405471216302666%3Fshowall%3Dtrue) avec BSEQ-sc) visent à ***déconvoluer*** les données bulk RNAseq pour estimer la proportion d'une population de cellules, **en apprenant à partir de l'abondance des transcrits scRNAseq**.

L'intérêt est que tout le monde aimerait profiter des avantages du scRNAseq au moyen d'une technique beaucoup moins chère et corroborée.

Le but du projet est donc d'***offrir une approche sans référence*** pour réaliser quelque chose qui est généralement fait en utilisant une approche basée sur des données de référence.

L'idée est simple et divisée en deux parties :
- **Apprendre** : L'utilisateur peut fournir les données brutes scRNAseq avec la liste des *barcode-clusters* des cellules. Le programme diviser les données brutes en fichiers FASTQ (un pour chaque cluster) et génère une bdd de comptage de k-mers (comme sur iMOKA).
- **Déconvoluer** : lu'tilisateur fournit des données brutes bluk RNAseq et un jeu de clusters de comptage de k-mers. Le programme prend chaque read du bulk, les divise en k-mers et, basé sur l'occurence de chaque k-mer dans les différents clusters, assigne le read à différents clusters.

Les avantages de cette méthode est que l'utilisateur peut prendre différents clusters de différentes études scRNAseq et les combiner. La sortie sera un fichier FASTQ, qui pourra être utilisé en tant qu'entrée par n'importe quel pipeline d'analyse. Les reads non assignés peuvent se montrer intéressants et révéler des types cellulaires non observés, de nouveaux transcrits, etc.


## Étape par étape

Le programme le plus utilisé pour l'analyse scRNAseq est CellRanger ([10X Genomics](https://www.10xgenomics.com/)). Il prend le fichier brut FASTQ, corrige les barcodes et aligne avec le programme STAR les cellules uniques à un génome de référence avec un comptage de gènes.

Par défaut, il regroupe aussi les cellules en clusters, mais ces clusters sont rarement utilsés dans les publications, le package R Seurat lui étant préféré ([exemple](https://bioinformaticsworkbook.org/dataAnalysis/RNA-Seq/Single_Cell_RNAseq/Chromium_Cell_Ranger.html#gsc.tab=0)), parce qu'il permet de mieux filtrer les données et la clusterisation est plus personnalisable.

Sans utiliser CellRanger, il serait plus difficile de réaliser la correction des barcodes (et puisque le programme fonctionne bien à ce niveau, il n'est pas utile de refaire ce travail). Les informations à propos des barcodes sont fournies dans un fichier BAM qui contient tous les reads (dont les non mappés).

## Ce qui a déjà été réalisé par Claudio

### Jeu de données
D'abord, il a fallut chercher un jeu de données. On a repris le jeu de donnée fournit par l'équipe qui a dévelopé [Music](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6342984/), et utilisé les données sur la souris, puisque cela ne nécessite aucun formulaire d'autorisation. Les données de Park et Beckermann ont été téléchargée.

### Programmes

#### CellRanger
Ensuite, CellRanger a été utilisé sur un échantillon pour l'utiliser en tant que test (puis sur tous les échantillons, donc ils sont déjà tous prêts), en utilisant le génome de référence de la souris fournit par 10X Genomics.

#### Seurat
Seurat n'a pas encore été utilisé, mais cela devrait être la prochaine étape notamment parce que le papier de Park merge les 7 échantillons en un seul en utilisant `cell_ranger aggregate` pour les combiner, et ensuite clusterise toutes les cellules ayant une bonne qualité (voir [supplementary materials](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6342984/#)). Du coup, il doit y avoir des changements à réaliser dans le cas où un utilisateur souhaite utiliser plusieurs fichiers BAM (parce que le merge n'utiliser que le comptage de gènes) et un seul fichier tsv barcode-cluster (on ne sait pas encore comment il prend en charge la collision des barcodes au sein de différents échantillons).

### Code

En ce qui concerne le code, tout a été écrit en Python au départ, mais comme c'était trop lent, tout a été réécrit en tant que partie de iMOKA.

Le code C++ (dans la branche de développement de iMOKA) prend le fichier BAM et le parse (réutilisation du code IRFINDER mais réimplémenté de zéro parce que IRFinder prenait des informations qui ne sont pas nécessaires ici).

#### Machine Learning

##### Apprentissage

Enfin, la partie apprentissage fonctionne très bien, mais il sera nécessaire de l'ajuster pour pour obtenir plusieurs fichiers BAM et pour prendre en main le problème des données scRNAseq mergées. KMC3 a été utilisé pour décomposer les fichiers FASTQ du cluster en comptage de k-mer, les trier et enfin créer la bdd binaire de comptage de k-mer.

L'apprentissage en C++ prend deux fois moins de temps qu'en Python.

##### Décomposition / déconvolution

La partie décomposition prend beaucoup plus de temps.

Un objet python a été créé afin de lire et faire des requêtes dans les fichiers de la bdd iMOKA

Le code C++ prend en entrée le fichier FASTQ (ou les fichiers en cas de paired-end, mais ce n'est pas encore bien implémenté, puisque le test actuel est single end), et prend un read à la fois, le divise en k-mer, ignore les k-mers qui commencent ou se terminent avec `max(3, k_len/3)` nucléotides A ou T (pour éviter les queues poly-A), et fait un comtpage de ces k-mers à partir de chaque cluster de la bdd des k-mers.

La normalisation de ces comptages est basée sur le nombre total de cellules dans le cluster original.

Une fois normalisés, les reads sont augmentés par un facteur A définit par l'utilisateur (défaut : 100), ce qui signifie qu'on assigne 100 instances du même read à différents clusters.

Chaque cluster $j$ va recevoir `round(A * Pj)` reads, où $P_j$ est le comptage normalisé trouvé avant, divisé par la somme de tous les comptages normalisés de tous les clusters.

Pour plus de clarté, regarder le slide [Miro](https://miro.com/welcomeonboard/NDhJR2NtNXZWQlBSZ1BYYTlQSzBqMmRQUFZjSXRXM3d2Nnhrc2tWMU1vWmJwd25MQmtScEhBWkxlU25JQ1VYOXwzMDc0NDU3MzY3NjM0OTc0Nzk4?invite_link_id=861101796954).

##### bdd

Les bdd iMOKA sont chargées dans le module de requête, par conséquent, tous les préfixes sont chargés en mémoire (pas les suffixes ; une solution extrême pour réduire me temps de traitement serait de de tout charger en mémoire pour éviter les lectures sur le disque).

Comme cela prenait du temps, un paramètre `--shift` a été ajouté pour ne prendre en considération que les k-mers qui se chevauchent d'une valeur shift (`--shift 1` pour prendre tous les k-mers), et s'il n'en trouve aucun, va prendre en considération. Les premiers tests ont montré que `--shift 7` diminue le temps d'éxecution de moitié, mais d'autres tests peuvent être nécessaires pour déterminer si c'est précis ou pas.

D'autres tests seront certainement nécessaires pour déterminer une différente valeur de k. On peuse que 31 (la valeur utilisée pour les tests réalisés) est trop grande et que 23 pourrait être suffisant.
