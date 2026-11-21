*Niveau L3 (bac +3) · Langage C · Durée indicative : 5 séances de 4 h + travail personnel*

## 1. Objectifs et prérequis

**Objectif final** : un analyseur de spectre qui capte un signal audio et affiche en temps réel l’énergie par bande de fréquence (des barres qui montent et descendent), d’abord sur PC dans un terminal, puis sur une carte STM32 Nucleo avec un écran dont le choix reste à faire. La version finale doit être **optimisée** : rapide, sobre en mémoire, sans calcul inutile.

**Compétences visées**

- Comprendre et programmer la transformée de Fourier rapide (FFT), le fenêtrage, les échelles logarithmiques.
- Séparer un code **portable** (le cœur de calcul) du code **dépendant de la plateforme** (audio, affichage).
- Optimiser de façon *mesurée* : virgule flottante simple précision, virgule fixe, tables précalculées, bibliothèque CMSIS-DSP.

**Prérequis** : C (pointeurs, structures, tableaux, types entiers de `stdint.h`), notions de nombres complexes, notions de trigonométrie. Aucune connaissance préalable en traitement du signal n’est supposée : la partie 0 donne les bases.

**Matériel et logiciels**

- PC Linux, macOS ou Windows (WSL possible, mais l’audio y est plus délicat) avec `gcc` et PortAudio.
- Une carte **Nucleo à cœur Cortex-M4F** (par exemple Nucleo-F411RE, Nucleo-F446RE ou Nucleo-G474RE) et STM32CubeIDE (ou une chaîne `arm-none-eabi-gcc` + OpenOCD).
- Un micro (électret avec préamplificateur, ou micro MEMS) pour la partie embarquée.
- Optionnel mais très utile : Python avec NumPy pour générer des signaux de référence et comparer vos résultats.

**Méthode** : chaque partie se termine par un critère de validation **testable**. Ne passez pas à la suite tant que celui-ci n’est pas satisfait. Les parties 1 à 4 se font sur PC, les parties 5 à 7 visent la carte.

## 2. Architecture logicielle imposée

La clé de ce TP est de ne jamais mélanger le calcul et les entrées-sorties. Le même code de calcul doit tourner sur PC et sur microcontrôleur.

```
spectre/
├── core/                  portable, aucun include système
│   ├── fft.h / fft.c
│   ├── spec.h / spec.c
│   └── tables.h           fenêtre, facteurs de rotation
├── pc/
│   ├── main_pc.c
│   ├── audio_pc.c         PortAudio
│   └── display_term.c     barres ANSI
├── mcu/
│   ├── main_mcu.c
│   ├── audio_adc_dma.c    ADC + timer + DMA
│   └── display_xxx.c      écran à choisir (ou UART au début)
├── tests/
│   └── test_core.c        signaux synthétiques
└── Makefile
```

Trois interfaces seulement relient les blocs :

- **Source audio** : fournit des trames de `SPEC_N` échantillons `int16_t`.
- **Cœur** : transforme une trame en `SPEC_BANDS` niveaux `uint8_t` (0 à 255).
- **Affichage** : reçoit ces niveaux et les dessine. Il ne sait rien de la FFT.

**Règles pour le dossier `core/`** : pas de `malloc`, pas de `printf`, pas de `double`, pas d’inclusion de PortAudio ou de HAL. Toute la mémoire est fournie par l’appelant.

## 3. Partie 0 — Théorie indispensable

*À lire **avant** de coder. Répondez aux questions de fin de section : elles reviendront dans le rapport.*

### 3.1 Échantillonnage et théorème de Nyquist-Shannon

Un signal audio est une tension qui varie continûment. L’ordinateur ou le microcontrôleur n’en voit que des **échantillons**, mesurés à intervalles réguliers de période Te = 1/fs, où fs est la fréquence d’échantillonnage.

Le théorème de Nyquist-Shannon dit qu’on ne peut représenter sans ambiguïté que les fréquences inférieures à **fs/2**. Au-delà, elles se **replient** (aliasing) : un sinus à 20 kHz échantillonné à fs = 32 kHz apparaît comme un sinus à |20 − 32| = 12 kHz. Il est alors impossible de le distinguer d’un vrai 12 kHz.

Conséquence pratique : sur une carte, un filtre **passe-bas analogique** (anti-repliement) doit précéder le convertisseur. Dans ce TP, on choisira **fs = 32 kHz**, soit une bande utile jusqu’à 16 kHz, largement suffisante pour un affichage musical.

*Questions* : 
(1) Quelle est la bande utile pour fs = 44,1 kHz ? 
(2) Un sinus à 30 kHz est échantillonné à 44,1 kHz : à quelle fréquence apparaît-il ?

### 3.2 De la transformée de Fourier discrète (TFD) à la FFT

Pour N échantillons x\[n\], la TFD calcule pour k = 0 … N−1 :

> X\[k\] = Σ (n = 0 … N−1) x\[n\] · e^(−j·2π·k·n / N)

**Interprétation** : X\[k\] mesure à quel point le signal « ressemble » à une sinusoïde complexe de fréquence **k·fs/N**. Son module |X\[k\]| donne l’amplitude de cette composante, son argument donne la phase. Pour un affichage de spectre, seul le module nous intéresse.

**Signal réel et symétrie** : pour un signal réel, X\[N−k\] est le conjugué de X\[k\]. Les valeurs utiles sont donc k = 0 … N/2 (de 0 Hz à fs/2). La seconde moitié est redondante. On exploitera cela à la partie 5.

**Coût** : la formule directe demande N multiplications complexes pour chacun des N coefficients, donc **N²**. Pour N = 1024, plus d’un million d’opérations par trame.

**L’idée de la FFT (Cooley-Tukey)** : séparer les échantillons d’indice pair et impair. Si E\[k\] est la TFD des échantillons pairs et O\[k\] celle des impairs (chacune sur N/2 points), alors :

> X\[k\] = E\[k\] + W^k · O\[k\]\
X\[k + N/2\] = E\[k\] − W^k · O\[k\]  avec W = e^(−j·2π/N)

Cette paire d’équations s’appelle un **papillon** : une multiplication complexe, une addition, une soustraction, pour deux sorties. En appliquant l’idée récursivement jusqu’à des TFD de taille 1, on obtient (N/2)·log2(N) papillons. Pour N = 1024 : **5120 papillons** contre 1 048 576 opérations, soit un facteur ≈ 200.

**Version itérative (celle qu’on programme)** :

1. Réordonner les échantillons par **inversion de bits** de l’indice (l’indice 1 = 0001 va en position 1000 = 8 pour N = 16). C’est l’ordre dans lequel les pairs/impairs se retrouvent après toutes les séparations.
2. Pour chaque étage de taille len = 2, 4, 8 … N, appliquer les papillons sur des blocs de len points, avec le facteur de rotation W\_len^j.

*Questions* : (1) Pourquoi N doit-il être une puissance de 2 pour cette version ? (2) Combien d’étages pour N = 1024 ? (3) Que vaut le bit-reverse de 3 sur 4 bits ?

### 3.3 Fenêtrage et fuite spectrale

La FFT suppose implicitement que la trame de N points se répète à l’infini. Si la trame ne contient pas un nombre entier de périodes du signal, la répétition crée une discontinuité aux bords, et l’énergie d’une raie « fuit » sur les bins voisins : c’est la **fuite spectrale**. Mathématiquement, tronquer le signal revient à le multiplier par une fenêtre rectangulaire, dont le spectre (un sinus cardinal) a de gros lobes secondaires.

On atténue ce défaut en multipliant la trame par une fenêtre qui tend vers 0 aux bords, par exemple la **fenêtre de Hann** :

> w\[n\] = 0,5 · (1 − cos(2π·n / (N−1)))

Compromis : lobes secondaires bien plus faibles, mais lobe principal plus large (la résolution se dégrade légèrement). De plus, la fenêtre diminue l’amplitude moyenne : pour Hann, le **gain cohérent** vaut 0,5.

*Question* : pourquoi calculer la fenêtre une fois pour toutes plutôt qu’à chaque trame ?

### 3.4 Résolution, échelle logarithmique, décibels

- **Résolution fréquentielle** : Δf = fs/N. Avec fs = 32 kHz et N = 1024, Δf = 31,25 Hz. Le bin k correspond à k·31,25 Hz.
- **Durée de la fenêtre** : N/fs = 32 ms. Plus N est grand, meilleure est la résolution en fréquence mais moins bonne la réactivité : c’est un compromis **temps-fréquence** incontournable.
- **Échelle logarithmique** : l’oreille perçoit les rapports de fréquence (une octave = un facteur 2), pas les différences. Une échelle linéaire gaspillerait presque tout l’écran sur les aigus. On répartit donc les bandes de façon logarithmique : bords f\_c = F\_MIN · (F\_MAX/F\_MIN)^(c/B) pour c = 0 … B.
- **Décibels** : niveau\_dB = 20·log10(|X| / référence), ou 10·log10(|X|²/référence²). Un sinus d’amplitude A fenêtré par Hann donne un pic de module **A·N/4** : c’est la normalisation à utiliser pour que 0 dB corresponde à un signal pleine échelle.
- **Dynamique** : un convertisseur 12 bits offre environ 72 dB, un 16 bits environ 96 dB. Afficher une plage de 60 à 80 dB est suffisant.

*Questions* : (1) Combien de bins se trouvent sous 100 Hz pour N = 1024, fs = 32 kHz ? Que cela implique-t-il pour les premières bandes logarithmiques ? (2) Pourquoi doubler N ne rend-il pas l’affichage « meilleur » dans tous les cas ?

### 3.5 Notions de temps réel : tampons, DMA, ping-pong

En temps réel, l’audio arrive en continu à cadence fixe : si votre traitement est trop lent, des échantillons sont perdus ou écrasés. Deux mécanismes :

- **Tampon circulaire** (anneau) : un tableau de taille puissance de 2 où l’on écrit à l’indice `w & (TAILLE-1)`. Le producteur (callback audio) écrit, le consommateur (boucle principale) lit les derniers échantillons, sans copie ni allocation.
- **DMA en double tampon (ping-pong)** : sur microcontrôleur, un contrôleur DMA recopie automatiquement les échantillons du convertisseur vers la RAM, sans mobiliser le processeur. Il lève une interruption à la **moitié** du tampon puis à sa **fin**. Pendant que le DMA remplit une moitié, le processeur traite l’autre.

**Contrainte fondamentale** : le traitement d’une moitié doit durer **moins que le temps de remplissage de l’autre moitié**. Avec une moitié de 512 échantillons à 32 kHz, cela fait 16 ms. C’est le « budget temps » que vous mesurerez à la partie 7.

**Recouvrement (overlap)** : pour une fenêtre de N = 1024 points renouvelée toutes les 512 nouvelles valeurs, on obtient un recouvrement de 50 %. L’affichage est alors rafraîchi toutes les 16 ms avec une fenêtre qui « glisse ».

*Questions* : (1) Pourquoi interdire `printf` et les mutex dans un callback temps réel ? (2) Que se passe-t-il si le traitement dure 20 ms au lieu de 16 ?

## 4. Partie 1 — Capture audio sur PC

**Objectif** : récupérer le flux du micro et le mettre à disposition du programme.

**Travail demandé**

1. Installer PortAudio (`sudo apt install portaudio19-dev` ou `brew install portaudio`).
2. Écrire un programme qui liste les périphériques (`Pa_GetDeviceCount`, `Pa_GetDeviceInfo`).
3. Ouvrir un flux d’entrée mono, `paFloat32`, fs = 32 000 Hz (ou 44 100 si le périphérique refuse), avec des blocs de 256 échantillons et un **callback**.
4. Dans le callback, copier les échantillons dans un tampon circulaire de 8192 floats. Écrire un index d’écriture de type `atomic_uint` (`<stdatomic.h>`).
5. Dans la boucle principale, toutes les 100 ms, lire les 1024 derniers échantillons, calculer leur **valeur efficace** (RMS = racine de la moyenne des carrés) et l’afficher.

**Validation** : RMS proche de 0 en silence, qui monte nettement quand vous parlez ou tapez dans les mains. Aucun `printf` dans le callback.

**Questions** : pourquoi la taille de l’anneau est-elle une puissance de 2 (que remplace-t-on alors ?). Que se passe-t-il à la lecture si l’index d’écriture dépasse 2³² ?

**Indices** : le callback reçoit un pointeur `const void *input` à convertir en `const float *`. Pour lire « les N derniers », partez de `w - N` en arithmétique non signée.

## 5. Partie 2 — La FFT

**Objectif** : programmer la FFT et la valider sur des signaux dont on connaît la réponse.

**Travail demandé**

1. Écrire `dft_naive` (double boucle, complexes en `double`) et la tester pour N = 256.
2. Écrire `fft` itérative radix-2 en place : inversion de bits puis papillons (cf. 3.2). Signature suggérée : `void fft(float complex *a, int n)`.
3. Comparer `fft` à `dft_naive` sur un signal aléatoire : l’écart maximal doit être négligeable (≈ 1e-9 en double, ≈ 1e-4 en float).
4. Mesurer le temps d’exécution des deux pour N = 256, 1024, 4096 (`clock_gettime`). Tracer ou tabuler : vérifier que le rapport suit la théorie (N² contre N·log N).
5. **Test du sinus** : générer x\[n\] = A·sin(2π·f·n/fs) avec fs = 32 000, N = 1024, f = 1000 Hz. Le pic doit se trouver au bin 1000·1024/32000 = **32** exactement.
6. Recommencer avec f = 1015 Hz (non multiple de Δf) : observer la fuite spectrale.
7. Appliquer la fenêtre de Hann et comparer. Vérifier que, pour f = 1000 Hz, |X\[32\]| ≈ A·N/4.

**Validation** : pic au bon bin, pic de module conforme à A·N/4, lobes secondaires nettement plus bas avec Hann pour le cas 1015 Hz.

**Questions** : pourquoi un sinus à fréquence exacte de bin ne produit-il pas de fuite sans fenêtre ? Où apparaît l’image « miroir » du pic dans le tableau de sortie ?

**Indices** : la boucle d’inversion de bits tient en cinq lignes (cherchez « bit-reversal permutation »). Le facteur de rotation de l’étage `len` est `cexpf(-2*I*M_PI/len)`, mis à jour par multiplication dans la boucle interne. Ne le recalculez pas avec `cos`/`sin` à chaque papillon : on le remplacera par une table à la partie 5.

## 6. Partie 3 — Du spectre à l’affichage terminal

**Objectif** : transformer les bins en barres lisibles, avec le même rendu que le premier programme de démonstration.

**Travail demandé**

1. **Bandes logarithmiques** : pour B colonnes (partez de 32 puis adaptez à la largeur du terminal), calculer les bords en Hz puis les bins : `lo = floor(f0·N/fs)`, `hi = ceil(f1·N/fs)`. Imposer `lo ≥ 1` (on ignore le continu) et `hi > lo`.
2. **Niveau d’une bande** : maximum des modules sur les bins de la bande (le maximum rend mieux les raies étroites que la moyenne).
3. **Conversion en dB** avec la normalisation A·N/4, puis mise à l’échelle linéaire de \[−80 dB ; 0 dB\] vers \[0 ; 1\].
4. **Lissage visuel** : attaque instantanée (si la nouvelle valeur est plus haute, on la prend), retombée lente (on retire une petite quantité à chaque image).
5. **Dessin** : séquences ANSI (`ESC[H` pour revenir en haut à gauche, `ESC[2J` pour effacer, `ESC[?25l` pour cacher le curseur), caractère plein `█` pour les barres. Construire toute l’image dans un tampon puis un **unique** `fwrite` pour éviter le scintillement.
6. Capter la taille du terminal avec `ioctl(TIOCGWINSZ)`.

**Validation** : jouez une tonalité pure depuis un générateur en ligne et balayez-la de 100 Hz à 10 kHz : la barre doit se déplacer de gauche à droite de façon régulière, à échelle logarithmique. Un bruit blanc doit donner des barres à peu près égales (en réalité légèrement croissantes vers la droite, car les bandes hautes contiennent plus de bins : à expliquer).

**Questions** : pourquoi les premières colonnes sont-elles souvent identiques entre elles ? Quelle est la conséquence d’un `DB_FLOOR` trop haut ? Trop bas ?

## 7. Partie 4 — Isoler le cœur portable

**Objectif** : obtenir un module `core/` testé, qu’on pourra compiler à l’identique pour le microcontrôleur.

**Travail demandé** : refactoriser votre code selon l’architecture de la section 2. Interface suggérée :

```c
// core/spec.h
#include <stdint.h>

#define SPEC_N      1024   // taille de la FFT
#define SPEC_BANDS  32     // nombre de barres
#define SPEC_FS     32000  // Hz

typedef struct {
    float    work[2 * SPEC_N];      // tampon de calcul complexe (re, im entrelacés)
    uint16_t band_lo[SPEC_BANDS];   // premier bin de chaque bande
    uint16_t band_hi[SPEC_BANDS];   // dernier bin + 1
    float    smooth[SPEC_BANDS];    // niveaux lissés
    uint8_t  level[SPEC_BANDS];     // sortie : 0..255
} spec_t;

void spec_init(spec_t *s);
void spec_process(spec_t *s, const int16_t *frame);  // frame : SPEC_N échantillons
const uint8_t *spec_levels(const spec_t *s);
```

Et deux interfaces pour les périphériques, que chaque plateforme implémente :

```c
// affichage
void display_init(void);
void display_draw_bars(const uint8_t *levels, int n);

// source audio : renvoie 1 quand une nouvelle trame est prête
int  audio_get_frame(int16_t *frame);
```

Écrire `tests/test_core.c` : génère des sinus de fréquences connues (100 Hz, 1 kHz, 5 kHz, 12 kHz), appelle `spec_process`, et vérifie (`assert`) que la bande la plus haute est celle attendue. Ajouter un test de silence (tout à zéro, aucune barre) et de signal pleine échelle (barre maximale).

**Validation** : `make test` passe ; `make pc` donne le même rendu qu’à la partie 3. Compiler `core/` avec `-Wall -Wextra -Wconversion -Wdouble-promotion` : zéro avertissement.

**Pourquoi des entrées `int16_t` et des sorties `uint8_t` ?** Ce sont exactement les formats naturels d’un ADC/I2S et d’un écran. Le code ne dépend plus du PC.

## 8. Partie 5 — Optimisation

**Principe directeur : mesurer avant d’optimiser.** Une optimisation non mesurée est une hypothèse. Pour chaque technique, relevez le temps (PC) ou les cycles (carte) **avant et après**, et consignez-les dans le tableau de la partie 7.

### 8.1 Les techniques à étudier

| # | Technique | Idée | Gain attendu |
| --- | --- | --- | --- |
| 1 | `float` au lieu de `double` | Le Cortex-M4F a une FPU simple précision matérielle. Le `double` est émulé en logiciel : 10 à 50 fois plus lent. Attention aux constantes (`1.0f`) et aux fonctions (`sinf`, `log10f`). | Très fort |
| 2 | Tables précalculées | Fenêtre de Hann et facteurs de rotation stockés en `const` (en Flash) au lieu de `cos`/`sin` dans la boucle. | Fort |
| 3 | FFT réelle | Une entrée réelle de N points se traite avec une FFT complexe de **N/2** points suivie d’une étape de recombinaison. | ≈ 2× |
| 4 | Pas de racine carrée | Travailler sur la puissance | X |
| 5 | Logarithme rapide | Approcher log2 par un comptage de zéros de tête (`__builtin_clz`, instruction CLZ) + petite table ou polynôme, puis multiplier par 3,0103 pour obtenir des dB. | Moyen |
| 6 | Virgule fixe Q15 | Entiers 16 bits représentant des réels dans \[−1 ; 1\[. Voir 8.2. | Variable, à mesurer |
| 7 | CMSIS-DSP | `arm_rfft_fast_f32`, `arm_rfft_q15` : FFT optimisées par ARM pour Cortex-M. | Référence à battre ou à adopter |
| 8 | Mémoire | Calcul en place, fenêtrage effectué pendant la copie des échantillons, aucune copie superflue. | Faible à moyen |
| 9 | Options de compilation | `-O2` ou `-Os`, `-mcpu=cortex-m4 -mfpu=fpv4-sp-d16 -mfloat-abi=hard`. Discuter `-ffast-math` (gain réel, mais ne respecte plus strictement IEEE). | Moyen |

### 8.2 La virgule fixe expliquée (notion délicate)

En Q15, un réel x de \[−1 ; 1\[ est stocké sous la forme de l’entier `round(x · 32768)`. Exemple : 0,5 → 16384.

- **Addition** : directe, mais attention au **débordement** : 0,75 + 0,75 dépasse la plage. On utilise alors une addition **saturée** (instruction `SSAT`/`QADD16` sur Cortex-M4) qui plafonne au lieu de s’enrouler.
- **Multiplication** : Q15 × Q15 donne un résultat en Q30 sur 32 bits. Il faut décaler de 15 bits à droite pour revenir en Q15 : `(int32_t)a * b >> 15`.
- **Dans la FFT** : chaque étage peut doubler l’amplitude (a + b). Pour éviter les débordements, on divise par 2 à chaque étage : après log2(N) étages, le résultat est divisé par N. Corollaire : la précision est limitée par le bruit de quantification, qui se relève face aux petits signaux. Les bandes faibles perdent quelques dB de dynamique.
- **Avantage** : sur Cortex-M4, les instructions DSP (`SMLAD`, `SMUAD`) font deux multiplications 16 bits par cycle ; la FFT Q15 est souvent la plus rapide et la plus économe en mémoire. Mais le M4F possède déjà une FPU : le flottant simple précision n’est pas pénalisant. **Il faut donc mesurer pour trancher**, ce qui est précisément l’exercice.

### 8.3 Travail demandé

1. Passer `core/` en `float`/`int16_t` purs (zéro `double`, vérifié par `-Wdouble-promotion`).
2. Générer la fenêtre et les facteurs de rotation par un petit programme (C ou Python) qui produit `tables.h` ; les inclure en `static const`.
3. Implémenter la FFT réelle via FFT complexe N/2 (cherchez « real FFT via complex FFT of half size » ; formule de recombinaison avec les facteurs W\_N^k).
4. Remplacer `sqrt` et `log10f` par la puissance et un log2 approché ; comparer l’erreur en dB sur une rampe d’amplitudes (elle doit rester inférieure à ≈ 0,5 dB, invisible sur l’écran).
5. Écrire une variante Q15 de la FFT (ou utiliser celle de CMSIS-DSP) et comparer son SNR à la version `float` sur vos sinus de test.
6. Documenter **pour chaque étape** le gain mesuré.

**Validation** : `make test` passe toujours (les tests de la partie 4 ne doivent pas changer), la version optimisée est au moins 3 fois plus rapide que votre première version float complexe, et l’écart de niveau entre version optimisée et version de référence reste inférieur à 1 dB sur les signaux de test.

**Questions** : (1) Pourquoi la FFT réelle n’apporte-t-elle pas exactement un facteur 2 ? (2) Que se passe-t-il en Q15 si un signal fort sature l’entrée ? (3) Dans quel cas la bibliothèque CMSIS-DSP ne vaut-elle pas le coup ?

## 9. Partie 6 — Portage sur STM32 Nucleo

**Objectif** : faire tourner `core/` sur la carte, avec une vraie source audio et un premier affichage de test.

### 9.1 Choix de la carte et de la source audio

- **Carte** : Cortex-M4F, au moins 64 Ko de RAM. Le cœur utilise environ 8 Ko (tampon complexe N = 1024 en float) plus quelques centaines d’octets : largement dans les moyens d’un F411RE, F446RE ou G474RE.
- **Option A (recommandée pour débuter)** : micro électret avec préamplificateur (module type MAX4466) sur une entrée de l’**ADC**. Le signal est centré sur Vcc/2 : il faudra **soustraire la valeur moyenne** (≈ 2048 en 12 bits) pour obtenir un signal centré sur zéro.
- **Option B (plus propre, plus complexe)** : micro MEMS I2S (type INMP441) lu par le périphérique **I2S** avec DMA. Les échantillons sont fournis sur 24 bits dans des mots de 32 bits : à décaler et convertir en `int16_t`.
- **Anti-repliement** : un petit filtre RC passe-bas (fc ≈ 10 à 15 kHz) avant l’ADC est fortement conseillé. Calculer R et C.

### 9.2 Chaîne d’acquisition (option A)

1. **Timer** configuré pour déclencher la conversion à exactement 32 kHz (sortie TRGO du timer connectée à l’ADC). L’échantillonnage doit être régulier : on ne le confie jamais à une boucle logicielle.
2. **ADC** en mode déclenché par le timer, conversions de 12 bits.
3. **DMA** en mode **circulaire** vers un tableau `uint16_t dma_buf[1024]` (deux moitiés de 512).
4. Callbacks HAL `HAL_ADC_ConvHalfCpltCallback` et `HAL_ADC_ConvCpltCallback` : ils se contentent de **lever un indicateur** `volatile` (`half_ready`, `full_ready`). Aucun calcul dans l’interruption.
5. Boucle principale : quand un indicateur est levé, reconstruire la trame de N = 1024 points, soustraire l’offset, convertir en `int16_t` avec fenêtrage au passage, appeler `spec_process`, puis `display_draw_bars`.

**Point d’attention (recouvrement)** : le tampon DMA contient toujours les 1024 derniers échantillons, mais dans un ordre qui dépend de la moitié qui vient d’être remplie. Quand la **première moitié** vient d’être complétée, les échantillons les plus anciens sont dans la seconde : la trame chronologique est donc *\[seconde moitié, première moitié\]*. Inversement à l’interruption de fin de tampon. *Question : pourquoi est-ce important pour la FFT ?*

### 9.3 L’écran : une décision différée, une architecture prête

L’écran n’est pas encore choisi : c’est précisément pour cela que seule la fonction `display_draw_bars(levels, n)` le connaît. Pistes à comparer dans une courte étude comparative (prix, bande passante, difficulté du pilote, lisibilité) :

| Type | Exemple | Remarques |
| --- | --- | --- |
| OLED monochrome I2C/SPI | SSD1306 128×64 | Simple et bon marché. 128 colonnes : 32 bandes de 4 px. L’I2C est lent : privilégier le SPI ou ne redessiner que ce qui change. |
| TFT couleur SPI | ST7789, ILI9341 | Rendu attrayant. Bande passante élevée : un redessin complet à 30 images/s est coûteux ; ne mettre à jour que la **différence** de hauteur de chaque barre, utiliser le DMA SPI. |
| Matrice de LED / ruban adressable | WS2812, MAX7219 | Très adaptée aux barres. Contrainte de timing stricte pour le WS2812. |

**Stratégie de test sans écran** : implémenter `display_uart.c`, qui envoie par l’UART de la Nucleo (accessible via l’USB du ST-LINK) le même rendu ANSI que `display_term.c`. Un terminal série (`picocom`, PuTTY) sur le PC permet alors de **valider toute la chaîne** avant d’avoir l’écran. Attention : l’envoi de plusieurs centaines d’octets par image à 115 200 bauds prend du temps (≈ 87 µs/octet) ; passer à 921 600 bauds et/ou n’envoyer que les barres qui ont changé.

**Contrainte d’architecture** : l’affichage ne doit **jamais** retarder l’acquisition. S’il est plus lent que les 16 ms d’une moitié de tampon, découpler : le traitement continue à chaque trame, l’affichage ne prend que le dernier résultat disponible (30 images/s suffisent).

### 9.4 Check-list de portage

- Horloge système et FPU activées (`SCB->CPACR`, normalement faite par le code de démarrage CubeIDE avec le bon réglage de l’option flottante).
- Taille de pile suffisante (`_Min_Stack_Size` dans le script d’édition de liens) ; ne placez pas `spec_t` sur la pile : variable globale ou statique.
- Retargeting de `printf` vers l’UART **uniquement pour le debug**, jamais dans le chemin temps réel.
- Même `tables.h` et même `core/` que sur le PC : aucune copie divergente.

**Validation** : avec le générateur de tonalité du PC placé devant le micro, la barre correspondante apparaît au bon endroit sur l’affichage test (UART). Pas de perte d’échantillons : une LED témoin s’allume si le traitement d’une trame n’est pas terminé quand la suivante arrive (détection de dépassement).

## 10. Partie 7 — Mesures, validation et rapport

### 10.1 Mesurer sur la carte

Le cœur Cortex-M4 contient un compteur de cycles matériel (DWT). Activation : `CoreDebug->DEMCR |= CoreDebug_DEMCR_TRCENA_Msk; DWT->CYCCNT = 0; DWT->CTRL |= DWT_CTRL_CYCCNTENA_Msk;`. Lire `DWT->CYCCNT` avant et après la fonction à mesurer. Conversion : temps = cycles / fréquence\_CPU.

Calculer ensuite : **charge CPU** = temps de traitement d’une trame / 16 ms.

### 10.2 Tableau de synthèse à compléter

| Version | Cycles FFT | Cycles par trame (total) | Charge CPU | RAM | Flash | Erreur vs référence |
| --- | --: | --: | --: | --: | --: | --: |
| V0 : double, FFT complexe, `cos`/`sin` |  |  |  |  |  |  |
| V1 : float |  |  |  |  |  |  |
| V2 : + tables |  |  |  |  |  |  |
| V3 : + FFT réelle |  |  |  |  |  |  |
| V4 : + puissance et log rapide |  |  |  |  |  |  |
| V5 : Q15 ou CMSIS-DSP |  |  |  |  |  |  |

### 10.3 Validation croisée PC / carte

Générer un signal synthétique (somme de trois sinus) stocké en Flash. Faire calculer les niveaux par la carte et par le PC avec le **même** `core/`, puis comparer les deux vecteurs (envoi par UART). Ils doivent être **identiques** pour une version flottante, et à ±1 niveau près pour la version Q15.

### 10.4 Rapport (6 à 8 pages)

1. Architecture et justification (schéma des blocs et des interfaces).
2. Théorie utile en vos propres mots : FFT, fenêtrage, dB (une demi-page chacun).
3. Résultats des validations des parties 2 et 3 (captures, courbes).
4. Tableau de la section 10.2, **commenté** : quelle optimisation a le meilleur rapport gain/effort ? Laquelle vous a surpris ?
5. Choix de l’écran et justification.
6. Limites connues et pistes d’amélioration.

## 11. Livrables et évaluation

- Dépôt Git avec historique de commits lisible (au moins un commit par partie).
- Le code complet, qui compile avec `make test`, `make pc` et le projet de la carte.
- Le rapport.
- Une démonstration sur la carte.

| Critère | Poids |
| --- | --: |
| Correction du cœur (FFT, bandes, dB) et tests | 25 % |
| Architecture portable et propreté du code | 20 % |
| Fonctionnement temps réel sur la carte (pas de perte, charge mesurée) | 20 % |
| Qualité de la démarche d’optimisation (mesures, analyse) | 20 % |
| Rapport et démonstration | 15 % |

## 12. Pour aller plus loin

- **Peak hold** : petite ligne qui retient le maximum de chaque barre avant de retomber.
- **Gain automatique** : adapter la plage dB au niveau moyen du signal.
- **Spectrogramme** : afficher l’évolution temporelle sur un écran couleur (défilement vertical).
- **Algorithme de Goertzel** : si seules quelques fréquences précises comptent (par exemple 8 bandes), calculer directement ces bandes sans FFT complète peut coûter moins cher. À comparer.
- **Filtres numériques (biquads)** avant l’analyse, pour appliquer un égaliseur et visualiser son effet.
- **Sortie audio** par DAC ou I2S pour entendre le signal filtré pendant que l’on voit son spectre.

## Annexe — Formules à garder sous la main

| Grandeur | Formule |
| --- | --- |
| Fréquence du bin k | f\_k = k · fs / N |
| Bin d’une fréquence f | k = f · N / fs |
| Résolution | Δf = fs / N |
| Durée de la fenêtre | T = N / fs |
| Fenêtre de Hann | w\[n\] = 0,5 · (1 − cos(2π·n / (N−1))) |
| Module d’un sinus d’amplitude A (Hann) | A · N / 4 |
| Niveau en dB (puissance P = re² + im²) | 10 · log10(P / (A\_ref · N/4)²) |
| Bords de bande c sur B bandes | F\_MIN · (F\_MAX / F\_MIN)^(c / B) |
| Nombre de papillons | (N / 2) · log2(N) |
| Budget temps d’une moitié de tampon DMA | (taille\_demi\_tampon) / fs |
