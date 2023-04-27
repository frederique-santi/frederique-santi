##### Simuler nb_sample fois une (micro-)forêt et (micro-)bosquets
# janvier-février 2023
# Auteurs :
# Rémy Margueritat <remy38m@gmail.com>
# Frederique Santi <frederique.santi@inrae.fr>
# Xabi Cazenave <Xabi.Cazenave@inrae.fr>

#### Définition des caractéristiques

# Trouver diviseurs  https://calculis.net/diviseur
# N doit être divisible par 16 (=taille des bosquets)
# N varie par multiples de 512, car 512=16*32, mais on peut diviser 
# la population en 2 encore : 256
n<-12 # varie : 0.5, 1, 2 ... 20
N<-n*512 # nombre total de plants
#  n  .5	 1	 2	 3	 4	 5	 6	 7	 8	 9	10	11	12	 13  14 15 	 16	 17	18	 19	 20
#lig  16	16	32	32	32	40	48	56	64	64	64	64	64	 64	 64 80	 64	 64	96	 64	 80
#col  16	32	32	48	64	64	64	64	64	72	80	88	96	104	112	96	128	136	96	152	128
nrowmat<-64 # nombre de lignes, déduire du tableau ci-dessus
ncolmat<-96 # nombre de colonnes, déduire du tableau ci-dessus
plantespece<-n*16 # nombre de plants par espèce
dp<-0.6 # distance entre plants en mètres
densite <- 1/(dp*dp) # densité de plants au mètre carré, plant au centre de dp*dp 
dm<-8 # distance minimale entre dominants en mètres, dm varie au cours du temps

# Déterminer le nombre d'espèces par catégorie de hauteur, 
# Dont 1/4 pour 1 = arbustes
# Puis calcul du nombre de plants par catégorie
especetotal<-32 # nombre d'espèces 8+6+6+6+6
espece1<-8 ; espece2<-6 ; espece3<-6 ; espece4<-6 ; espece5<-6
categorie1<-plantespece*espece1 ; 
categorie2<-plantespece*espece2 ; 
categorie3<-plantespece*espece3 ; 
categorie4<-plantespece*espece4 ; 
categorie5<-plantespece*espece5

# Attribuer des numéros de plants dans chaque espèce/catégorie
# Les numéros d'espèces sont multipliés par 1000 
# Le nombre de plants pour chaque catégorie est inférieur à 1000
plant<-NULL
for (i in 1:especetotal) plant<-c(plant,c((1:plantespece)+i*1000))

# catégorie de hauteur, quand plant est en ordre croissant
categorie<-c(rep(1,categorie1),rep(2,categorie2),rep(3,categorie3),
             rep(4,categorie4),rep(5,categorie5))

# numérotation dans le sens de lecture du processus itératif à venir
numero<-matrix(1:N,nrow=nrowmat,byrow=TRUE)
numero_encolonne<-matrix(numero)

# coordonnées xx (lignes) et yy (colonnes) en mètres
xx<-rep(seq(0,dp*(nrowmat-1),by=dp),ncolmat) 
yy<-rep(seq(0,dp*(ncolmat-1),by=dp),each=nrowmat)

# statut à 1 = dominé par défaut
statut<-rep(1,N)

#### Gérération de 100 micro-forêts et micro-bosquets équivalents
# chaque ensemble de dominants sera placé dans une liste
result_list_d <- list()
result_list_e <- list()
nb_sample <- 1000
result_de <- matrix(data = NA, nrow = nb_sample, ncol = 7, byrow = FALSE)
       dimnames(result_de)[[2]] = c("rep","nb_dom_d","nb_esp_d","nb_cat_d",
                    "nb_dom_e","nb_esp_e","nb_cat_e")

for (j in 1:nb_sample){
    set.seed(j) # on initialise le tirage pour avoir un résultat reproductible

# Génération d'une micro-forêt
a<-sample(plant)
b<-matrix(a,nrowmat,ncolmat)

### Décrire la micro-foret avec dataframe = une ligne un plant
# avec numero, plant, espece, coordonnées xx et yy, catégorie,
# statut (1 = dominé, 1000 = dominant pour que ce soit visible), 
# classe qui décrit la compétition inter-espèces d'une même catégorie

# espece
espece<-floor(a/1000)

# Fusion avec numero espece plant xx yy
d<-cbind(numero_encolonne,espece,a,xx,yy) # fusion en colonne
colnames(d)<-c("numero","espece","plant","xx","yy")
d<-as.data.frame(d)

# Fusion avec la catégorie, la classe aléatoire intra catégorie et le statut
d<-d[order(d$plant),] # tri en ordre croissant des numéros de plants
classe<-matrix(c(rep(sample(1*1000+1:categorie1,replace=FALSE)),
                 rep(sample(2*1000+1:categorie2,replace=FALSE)),
                 rep(sample(3*1000+1:categorie3,replace=FALSE)),
                 rep(sample(4*1000+1:categorie4,replace=FALSE)),
                 rep(sample(5*1000+1:categorie5,replace=FALSE))))
# Rajouter un diamètre de houppier
# l'ensemble des arbres et arbustes suivraient une loi normale
# l'écart type dépend de la composition entre espèces, prendre sd=1/4*mean
meanhouppier<-5
houppier<-sort(rnorm(n*512, mean = meanhouppier, sd = 0.25*meanhouppier))
d<-cbind(d,categorie,classe,houppier,statut)

### Trouver les dominants 1000
for (i in c(1:N)) {
  current_value <- d$classe[i]
  is_neighbor <- sqrt((d$xx - d$xx[i])^2 + (d$yy - d$yy[i])^2) < dm
  is_neighbor[i] <- FALSE
  neighbors <- d$classe[is_neighbor]
  d$statut[i] <- ifelse(current_value > max(neighbors),
                        1000,d$statut[i])
}
d_dom <- d[d$statut==1000,]

# Visualisation intermédiaire de la classe sous la forme de b
# b_classe <- matrix(d$statut, nrow = nrowmat, ncol = ncolmat, byrow=TRUE) 

### Eclater la foret en bosquets et statut= 1 dominés
# on rajoute des coordonnées éclatées xxx et yyy (exemple + 10 m)
e<-d
e<-e[order(e$numero),]  # données rangées col par col /e$xx et e$yy
# mais on veut créer e$xxx d'abord
e<-e[order(e$yy),]  # données rangées ligne par ligne
e$xxx<- rep(c(0,dp,2*dp,3*dp), times=ncolmat*nrowmat/4) +
  rep(c(0,0,0,0,rep(10*(1:(nrowmat/4-1)),each=4)), times=ncolmat)
e<-e[order(e$xx),]  # données rangées col par col
e$yyy<-rep(c(0,dp,2*dp,3*dp), times=nrowmat*ncolmat/4) +
          rep(c(0,0,0,0,rep(10*(1:(ncolmat/4-1)),each=4)), times=nrowmat)

# remettre le statut à 1 pour tous
e$statut <- 1

### Trouver les dominants 1000 des microbosquets e
for (i in c(1:N)) {
  current_value <- e$classe[i]
  is_neighbor <- sqrt((e$xxx - e$xxx[i])^2 + (e$yyy - e$yyy[i])^2) < dm
  is_neighbor[i] <- FALSE
  neighbors <- e$classe[is_neighbor]
  e$statut[i] <- ifelse(current_value > max(neighbors),
                        1000,e$statut[i])
}
e_dom <- e[e$statut==1000,]

# Visualisation intermédiaire de la classe sous la forme de b
# e_classe <- matrix(e$statut, nrow = nrowmat, ncol = ncolmat, byrow=TRUE) 

### remplir les list et matrix de résultats
result_list_d[[j]] <- d_dom
result_list_e[[j]] <- e_dom
# Compter le nombre d'arbres, d'especes, de categories de d_dom et e_dom
result_de[j,] <-c(j, dim(d_dom)[[1]],
                  length(unique(d_dom$espece)),
                  length(unique(d_dom$categorie)),
                  dim(e_dom)[[1]],
                  length(unique(e_dom$espece)),
                  length(unique(e_dom$categorie)))

} # terminer la boucle j

# sortir les fichiers résultats en faisant varier dm de 3 à 8
# options(width = 90)         # pour sample = 1000
# options(max.print = 300000) # pour sample = 1000
# sink("result_de_dm8_sample1000.txt") ; result_de ; sink()
# sink("result_de_dm8_data.txt") ; result_list_d ; result_list_e ; sink()

# sortir les fichiers résultats dm=8 donc avec canopée "finale"
# tailles de pops variables : n varie : 0.5, 1, 2 ... 20 
sink("result_de_n12_sample1000.txt") ; result_de ; sink()

#### Créer des bosquets avec 4 plants pour 3 espèces arborées différentes 
# et 1 plant par espèce d'arbuste (catégorie1), toutes différentes
# total 4 4 4 1 1 1 1 : 16 plants et 7 espèces par bosquet
# une des espèces arborée dominera et les autres seront éliminées naturellement
# ou par éclaircie volontaire
# la disposition peut être aléatoire, mais pour favoriser la sélection 
# intra-espèce, privilégier une disposition par 2 pour les espèces arborées 
# placer les arbustes dans les angles pour qu'ils persistent
# exemple :
# a1	A1	A2	a2
# A3	A1	A2	A1
# A3	A2	A3	A1
# a3	A2	A3	a4

### Les noms d'espèces en bonne quantité
#  nombre de bosquets : n*512/16=n*32
n<-2
bosquet4<-c(rep(9:32,times=8)) # avant sample

### Echantillonner 2300 fois les 3 espèces d'arbres des micro-bosquets
# cela permet d'éliminer ensuite (par excel) les cas où PB3>0, puis d'obtenir
# 1000 runs avec PB=0, par souci de comparaison avec le cas précédent
nb_sample<-2600 
# compo_bosquet4 <- list()
result_f <- matrix(data = NA, nrow = nb_sample, ncol = 5, byrow = TRUE)

for (l in 1:nb_sample){
  set.seed(l) # on initialise le tirage pour avoir un résultat reproductible
  
## le dominant de chaque bosquet et le nombre d'espèces différentes dominantes
# Premier sample
bosquet4<-sample(bosquet4,replace=FALSE)

## faire le test pour chaque bloc de 3 successivement
# mettre à part les blocs où une espèce est répétée 2 fois et resample
# traiter les trois cas de répartition entre catégories 
# de type 333, 345, 334 ou 433 (avec exemple des catégories 3, 4, 5)

f_dom<-NULL ;PB<-NULL ;  PB2<-NULL ;  PB3<-NULL ; OK<-NULL

for (k in seq(from=1, to=length(bosquet4), by=3)) {
  bloc<-sort(c(bosquet4[k],bosquet4[k+1],bosquet4[k+2]))
  # tester si une valeur est répétée, si oui mettre dans fichier "PB"
  if(bloc[1]-bloc[2]==0||bloc[2]-bloc[3]==0) PB<-c(PB,bloc) else {
    # si non mettre le bloc dans fichier "OK"
    # traiter cas 333 : prendre le premier avant sort puisque aléatoire
    if(bloc[3]-bloc[1]<6) { f_dom<-c(f_dom,bosquet4[k]); OK<-c(OK,bloc)} else { 
        # traiter cas 345 : prendre le max=3eme
        if(bloc[2]-bloc[1]>5 && bloc[3]-bloc[2]>5) { f_dom<-c(f_dom,bloc[3]) 
        OK<-c(OK,bloc)} else {           
            # traiter cas 334 : prendre le max=3eme
            if(bloc[2]-bloc[1]<6) {f_dom<-c(f_dom,bloc[3]);OK<-c(OK,bloc)} else{
              # traiter cas 344 : prendre le premier du sample 44
              f_dom<-c(f_dom,sample(c(bloc[2],bloc[3]))[1]) ; OK<-c(OK,bloc)
              }}}}}
              # fin 1ere boucle de test k

## resample puis tout refaire pour PB
if (length(PB)>0) { 
PB<-sample(PB)
for (k in seq(from=1, to=length(PB), by=3)) {
  bloc<-sort(c(PB[k],PB[k+1],PB[k+2]))
  # tester si une valeur est répétée, si oui mettre dans fichier "PB2"
  if(bloc[1]-bloc[2]==0||bloc[2]-bloc[3]==0) PB2<-c(PB2,bloc) else {
    # si non mettre le bloc dans fichier "OK"
    # traiter cas 333 : prendre le premier avant sort puisque aléatoire
    if(bloc[3]-bloc[1]<6) { f_dom<-c(f_dom,PB[k]); OK<-c(OK,bloc)} else { 
      # traiter cas 345 : prendre le max=3eme
      if(bloc[2]-bloc[1]>5 && bloc[3]-bloc[2]>5) { f_dom<-c(f_dom,bloc[3]) 
      OK<-c(OK,bloc)} else {           
        # traiter cas 334 : prendre le max=3eme
        if(bloc[2]-bloc[1]<6) {f_dom<-c(f_dom,bloc[3]);OK<-c(OK,bloc)} else{
          # traiter cas 344 : prendre le premier du sample 44
          f_dom<-c(f_dom,sample(c(bloc[2],bloc[3]))[1]) ; OK<-c(OK,bloc)
          }}}}}}
          # fin 2eme boucle de test k

## resample puis tout refaire pour PB2
if (length(PB2)>0) {
PB2<-sample(PB2)
for (k in seq(from=1, to=length(PB2), by=3)) {
  bloc<-sort(c(PB2[k],PB2[k+1],PB2[k+2]))
  # tester si une valeur est répétée, si oui mettre dans fichier "PB3"
  if(bloc[1]-bloc[2]==0||bloc[2]-bloc[3]==0) PB3<-c(PB3,bloc) else {
    # si non mettre le bloc dans fichier "OK"
    # traiter cas 333 : prendre le premier avant sort puisque aléatoire
    if(bloc[3]-bloc[1]<6) { f_dom<-c(f_dom,PB2[k]); OK<-c(OK,bloc)} else { 
      # traiter cas 345 : prendre le max=3eme
      if(bloc[2]-bloc[1]>5 && bloc[3]-bloc[2]>5) { f_dom<-c(f_dom,bloc[3]) 
      OK<-c(OK,bloc)} else {           
        # traiter cas 334 : prendre le max=3eme
        if(bloc[2]-bloc[1]<6) {f_dom<-c(f_dom,bloc[3]);OK<-c(OK,bloc)} else{
          # traiter cas 344 : prendre le premier du sample 44
          f_dom<-c(f_dom,sample(c(bloc[2],bloc[3]))[1]) ; OK<-c(OK,bloc)
          }}}}}}
          # fin 3eme boucle de test k

## si PB3 existe, on ne resample pas
      if (length(PB3)>0) { 
        for (k in seq(from=1, to=length(PB3), by=3))   {
          bloc<-sort(c(PB3[k],PB3[k+1],PB3[k+2]))
        # traiter cas 3=3=3 : prendre le premier avant sort puisque aléatoire
          if(bloc[3]-bloc[1]<6) { f_dom<-c(f_dom,PB3[k]); OK<-c(OK,bloc)} else { 
          # traiter cas 3=34 : prendre le max=3eme
          if(bloc[2]-bloc[1]<6) {f_dom<-c(f_dom,bloc[3]);OK<-c(OK,bloc)} else{
                # traiter cas 34=4 : prendre le premier du sample 44
                f_dom<-c(f_dom,sample(c(bloc[2],bloc[3]))[1]) ; OK<-c(OK,bloc)
          }}}}
          # fin 4eme boucle de test k

# test intermédiaire avant de faire la boucle l
# c(length(f_dom),length(OK)/3,length(PB)/3,length(PB2)/3,length(PB3)/3)

### remplir les fichiers résultats pour chaque l
# mettre liste des essences dominantes et les length des PB
result_f[l,] <- c(l,length(unique(f_dom)),
                  length(PB),length(PB2),length(PB3))

# bosquet4 dans la liste de sortie, sous format 3 col, trié par ligne
# if (length(PB3)>0) {
# compo_bosquet4[[l]]<- matrix(data=OK, nrow = n*32, byrow = TRUE)} else {
  # seuls ceux qui n'ont pas de PB3 sont sortis
  # compo_bosquet4[[l]]<- NULL
  # }

} # fin de la boucle l

### sortir le fichier résultat
colnames(result_f)<-c("numero","esp_f",
                      "nb_PB_test1","nb_PB_test2","nb_PB_test3")
# nécessaire d'augmenter la taille imprimable
# car message [ reached getOption("max.print") -- omitted 200 rows ]
options(max.print = 600000)
sink("result_f.txt") ; result_f ; sink()
# 1016 résultats avec nb_PB_test3=0
# sink("compo_bosquet4.txt") ; compo_bosquet4 ; sink()

### créer les 4 col d'arbustes, si possible différents
compo_bosquet1 <- list()

## nombre de combinaisons différentes classées de 4 éléments dans 8 = 70
# on n'utilise pas arbuste<-c(rep(1:8,times=32)) puis des sample
dim(combn(8, 4)) # ce n'est pas en liste, par défaut simplify = TRUE
# [1]  4 70

## choisir 6 à exclure de 70 pour avoir le même nombre d'arbres
# 1234 5678 1256 1278 3456 3478
# Tous les autres seront uniques et le total = 64 OK
compo_bosquet1<-combn(8, 4, simplify = FALSE)  # retourne une liste
compo_bosquet1[[70]]<- NULL
compo_bosquet1[[61]]<- NULL
compo_bosquet1[[56]]<- NULL
compo_bosquet1[[15]]<- NULL
compo_bosquet1[[10]]<- NULL
compo_bosquet1[[1]]<- NULL

#### Visualiser une pop et des bosquets équivalents par la canopée
# en ordonnant le plus gros est au-dessus
d<-d[order(d$classe),]
e<-e[order(e$classe),]

# mettre d en bas à droite de e en créant xxx et yyy pour d
max(e$xxx) ; max(e$yyy) # 71.8 pour N=1040, dp=0.6
d$xxx<-d$xx + max(e$xxx)+20
d$yyy<-d$yy 
max(d$xxx) ; max(d$yyy) # 110.4 et 18.6
echelle_def <-c(NA,NA,NA,NA,NA,NA,NA,8,NA,max(d$xxx),max(d$xxx)) # sert à caler l'échelle
de<-rbind(d,e,echelle_def)

# compiler e_dom à d (en bas à droite comme témoin) et d_dom (en haut à droite)
d_dom$xxx<-d_dom$xx + max(e$xxx)+20
d_dom$yyy<-d_dom$yy+(max(e$yyy)-max(d$yy))
de_dom<-rbind(e_dom,d,d_dom,echelle_def)

# fabriquer un f_dom avec PB3=0
f_domOK<-f_dom # permet de revenir en arrière
# sample 4 plants par espèce dominante de f_domOK et prendre le max de classe
# et le houppier correspondant de d
# prendre 4 plants successifs dans chaque espèce et choisir la classe max de d
# puis créer le fichier complet f_domt
# en collant d[,1:8] et les statut (1000 partout) et localisations de e_dom
# intermédiaire pour que les 3 dernières col soient déjà OK = arbres mélangés
f_dom<-sort(f_dom)
f_domt<-e_dom[order(sample(e_dom$numero, replace = FALSE)),] 

f_domclasse<-NULL
for(i in 1:64) {
  if (i==1 || f_dom[i]!=f_dom[i-1]) { k=0
    f_domclasse[i]<-max(d[d$plant==f_dom[i]*1000+1,7],
                        d[d$plant==f_dom[i]*1000+2,7],
                        d[d$plant==f_dom[i]*1000+3,7],
                        d[d$plant==f_dom[i]*1000+4,7])
     } else { k=k+4
    f_domclasse[i]<-max(d[d$plant==f_dom[i]*1000+1+k,7],
                        d[d$plant==f_dom[i]*1000+2+k,7],
                        d[d$plant==f_dom[i]*1000+3+k,7],
                        d[d$plant==f_dom[i]*1000+4+k,7]) } 
  f_domt[i,1:8]<-d[d$classe==f_domclasse[i],1:8] 
              }

# compiler f_domt à d (en bas à droite comme témoin) et d_dom (en haut à droite)
df_dom<-rbind(f_domt,d,d_dom,echelle_def)

# conserver les valeurs des exemples et export
exemple<-list()
d<-d[order(d$numero),] ; exemple$d<-d
e<-e[order(e$numero),] ; exemple$e<-e
d_dom<-d_dom[order(d_dom$numero),] ; exemple$d_dom<-d_dom
e_dom<-e_dom[order(e_dom$numero),] ; exemple$e_dom<-e_dom
exemple$bosquet4<-bosquet4 ; exemple$f_domOK<-f_domOK
exemple$f_domclasse<-f_domclasse ; exemple$f_domt<-f_domt
exemple$nb_espèces_def_dom<-c(length(unique(d_dom$espece)),
                 length(unique(e_dom$espece)),length(unique(f_domOK)))
options(width = 90) # pour que toutes les colonnes soient ensemble
sink("def_exemple1") ; exemple ; sink()

# plot - install ggplot2 dans onglet Packages
# zoom, modifier taille, clic droit save image as
print(ggplot(df_dom,aes(x=xxx,y=yyy,color=espece,size=houppier)) 
      + geom_point() + scale_color_gradientn(colours=rainbow(6)))
print(ggplot(de_dom,aes(x=xxx,y=yyy,color=espece,size=houppier)) 
      + geom_point() + scale_color_gradientn(colours=rainbow(6)))
