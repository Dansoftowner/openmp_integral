# 2D Gauss-integrál numerikus közelítése

Ez a projekt a kétdimenziós Gauss-függvény numerikus integrálását vizsgálja C nyelven,
szekvenciális és OpenMP-vel párhuzamosított megoldásokkal. A programok célja a

```text
f(x, y) = e^(-(x^2 + y^2))
```

függvény integráljának közelítése különböző módszerekkel, valamint a szekvenciális és
párhuzamos futási idők összehasonlítása.

A teljes síkon vett ismert eredmény:

```text
integral integral e^(-(x^2 + y^2)) dx dy = pi
```

A projektben szereplő programok ezt az összefüggést használják referenciaértékként,
illetve körlapra korlátozott integrál esetén az analitikus részeredményt számolják ki.

## Projekt felépítése

```text
.
|-- gauss2d_integral_seq.c
|-- gauss2d_integral_par.c
|-- gauss2d_circle_integral_seq.c
|-- gauss2d_circle_integral_par.c
|-- gauss2d_montecarlo_integral_seq.c
|-- gauss2d_montecarlo_integral_par.c
`-- README.md
```

## Programok rövid leírása

| Fájl | Típus | Módszer | Párhuzamosítás |
| --- | --- | --- | --- |
| `gauss2d_integral_seq.c` | szekvenciális | determinisztikus rácsos integrálás transzformált tartományon | nincs |
| `gauss2d_integral_par.c` | párhuzamos | determinisztikus rácsos integrálás transzformált tartományon | OpenMP |
| `gauss2d_circle_integral_seq.c` | szekvenciális | polárkoordinátás integrálás körlapon | nincs |
| `gauss2d_circle_integral_par.c` | párhuzamos | polárkoordinátás integrálás körlapon | OpenMP |
| `gauss2d_montecarlo_integral_seq.c` | szekvenciális | Monte Carlo mintavételezés transzformált tartományon | nincs |
| `gauss2d_montecarlo_integral_par.c` | párhuzamos | Monte Carlo mintavételezés transzformált tartományon | OpenMP |

## Matematikai háttér

### 1. Teljes síkbeli Gauss-integrál

A kétdimenziós Gauss-integrál:

```text
I = integral_{-inf}^{inf} integral_{-inf}^{inf} e^(-(x^2 + y^2)) dx dy
```

szétválasztható két egydimenziós Gauss-integrál szorzatára:

```text
I = (integral_{-inf}^{inf} e^(-x^2) dx)^2 = sqrt(pi)^2 = pi
```

Ezért a programok a kapott numerikus eredményt `pi` értékéhez hasonlítják.

### 2. Végtelen tartomány transzformálása

A teljes síkbeli integrál numerikus rácsozása közvetlenül nehéz, mert a tartomány
végtelen. A projekt a következő helyettesítést használja:

```text
x = tan(u)
y = tan(v)
u, v eleme [-pi/2, pi/2]
```

A transzformáció Jacobi-determinánsa:

```text
dx/du = 1 / cos^2(u)
dy/dv = 1 / cos^2(v)
```

Ez alapján a transzformált integrandus:

```text
g(u, v) = e^(-(tan^2(u) + tan^2(v))) / (cos^2(u) cos^2(v))
```

A `gauss2d_integral_*` és `gauss2d_montecarlo_integral_*` programok ezt a
transzformált függvényt integrálják a véges `[-pi/2, pi/2] x [-pi/2, pi/2]`
tartományon.

### 3. Körlapon vett integrál polárkoordinátákban

A körlapos változat az `R` sugarú kör területén integrál:

```text
I_R = integral_0^R integral_0^(2pi) e^(-r^2) r dphi dr
```

Itt az `r` szorzó a polárkoordinátás Jacobi-tényező. Az analitikus eredmény:

```text
I_R = pi * (1 - e^(-R^2))
```

A `gauss2d_circle_integral_*` programok ezt az értéket is kiírják, és ehhez mérik a
numerikus hibát.

## Numerikus módszerek

### Középpontos téglalap-szabály

A determinisztikus integrálos programok középpontos téglalap-szabályt használnak.
A tartományt egyenlő részekre osztják, és minden cellában a cella középpontjában
értékelik ki az integrandust.

A transzformált teljes síkbeli integrál esetén:

```text
dx = (b - a) / N
dy = (d - c) / M
```

ahol:

```text
a = c = -pi/2
b = d =  pi/2
```

Az összegzés:

```text
sum += g(a + (i + 0.5) * dx, c + (j + 0.5) * dy)
eredmeny = sum * dx * dy
```

A körlapos változatban:

```text
dr = R / N_r
dphi = 2pi / N_phi
```

Az összegzés:

```text
sum += e^(-r_mid^2) * r_mid
eredmeny = sum * dr * dphi
```

### Monte Carlo közelítés

A Monte Carlo programok véletlen pontokat generálnak a transzformált
`[-pi/2, pi/2] x [-pi/2, pi/2]` tartományon. Az integrál becslése:

```text
eredmeny = terulet * atlag(g(u, v))
```

ahol a transzformált tartomány területe:

```text
terulet = pi * pi
```

A Monte Carlo módszer előnye, hogy egyszerű és jól párhuzamosítható. Hátránya, hogy
a pontosság statisztikus jellegű: az eredmény futásonként változhat, és a hiba
általában lassabban csökken, mint a determinisztikus rácsos módszer esetében.

## Párhuzamosítás

A párhuzamos programok OpenMP-t használnak.

### Rácsos integrálás

A determinisztikus párhuzamos programokban a két egymásba ágyazott ciklus
összevonva kerül párhuzamosításra:

```c
#pragma omp parallel for collapse(2) reduction(+:sum)
```

Ennek jelentése:

- `parallel for`: a ciklus iterációit több szál között osztja szét;
- `collapse(2)`: a kétdimenziós ciklust egy nagy iterációs térré alakítja;
- `reduction(+:sum)`: minden szál saját részösszeget kap, majd ezek biztonságosan
  összegződnek.

### Monte Carlo

A párhuzamos Monte Carlo programban:

```c
#pragma omp parallel reduction(+:sum)
```

blokk indul, ezen belül pedig:

```c
#pragma omp for
```

osztja szét a mintákat a szálak között.

A szekvenciális Monte Carlo program `rand()`-ot használ, a párhuzamos változat pedig
`rand_r()`-t, mert a `rand()` globális állapotú, és párhuzamos használatnál nem jó
választás.

## Követelmények

A fordításhoz olyan C fordító kell, amely támogatja:

- a C standard matematikai könyvtárat;
- `clock_gettime` és `CLOCK_MONOTONIC` használatát;
- OpenMP-t a párhuzamos fájlokhoz.

Linuxon vagy WSL-ben tipikusan elég a GCC:

```bash
gcc --version
```

Windows alatt ajánlott WSL, MSYS2/MinGW vagy más olyan környezet, ahol az OpenMP és
a POSIX idő API elérhető.

## Fordítás

### Szekvenciális programok

```bash
gcc -O2 gauss2d_integral_seq.c -lm -o gauss2d_integral_seq
gcc -O2 gauss2d_circle_integral_seq.c -lm -o gauss2d_circle_integral_seq
gcc -O2 gauss2d_montecarlo_integral_seq.c -lm -o gauss2d_montecarlo_integral_seq
```

### Párhuzamos programok

```bash
gcc -O2 -fopenmp gauss2d_integral_par.c -lm -o gauss2d_integral_par
gcc -O2 -fopenmp gauss2d_circle_integral_par.c -lm -o gauss2d_circle_integral_par
gcc -O2 -fopenmp gauss2d_montecarlo_integral_par.c -lm -o gauss2d_montecarlo_integral_par
```

Megjegyzés: egyes rendszereken a `-lm` kapcsolót a parancs végén kell hagyni, mert
a matematikai könyvtár linkelése sorrendfüggő lehet.

## Futtatás

### Teljes síkbeli integrál, rácsos módszer

Szekvenciális:

```bash
./gauss2d_integral_seq <N> <M>
```

Párhuzamos:

```bash
./gauss2d_integral_par <N> <M>
```

Példa:

```bash
./gauss2d_integral_seq 1000 1000
./gauss2d_integral_par 1000 1000
```

Paraméterek:

- `N`: felosztások száma az első dimenzióban;
- `M`: felosztások száma a második dimenzióban.

### Körlapon vett integrál

Szekvenciális:

```bash
./gauss2d_circle_integral_seq <R_sugar> <N_r_felosztas> <N_phi_felosztas>
```

Párhuzamos:

```bash
./gauss2d_circle_integral_par <R_sugar> <N_r_felosztas> <N_phi_felosztas>
```

Példa:

```bash
./gauss2d_circle_integral_seq 5 1000 1000
./gauss2d_circle_integral_par 5 1000 1000
```

Paraméterek:

- `R_sugar`: a kör sugara;
- `N_r_felosztas`: sugárirányú felosztások száma;
- `N_phi_felosztas`: szögirányú felosztások száma.

### Monte Carlo integrál

Szekvenciális:

```bash
./gauss2d_montecarlo_integral_seq <N_mintak_szama>
```

Párhuzamos:

```bash
./gauss2d_montecarlo_integral_par <N_mintak_szama>
```

Példa:

```bash
./gauss2d_montecarlo_integral_seq 10000000
./gauss2d_montecarlo_integral_par 10000000
```

Paraméterek:

- `N_mintak_szama`: a véletlen mintapontok darabszáma.

## OpenMP szálak számának beállítása

Az OpenMP programok futásánál a szálak száma az `OMP_NUM_THREADS` környezeti
változóval állítható.

Linux, WSL vagy Git Bash:

```bash
export OMP_NUM_THREADS=4
./gauss2d_integral_par 2000 2000
```

PowerShell:

```powershell
$env:OMP_NUM_THREADS = "4"
.\gauss2d_integral_par.exe 2000 2000
```

## Kimenet értelmezése

A programok tipikusan a következő adatokat írják ki:

```text
Eredmény: ...
Hiba: ...
Idő: ...
```

A körlapos programok ezen felül kiírják:

```text
Analitikus: ...
```

Jelentés:

- `Eredmény`: a numerikusan számolt integrál;
- `Analitikus`: a zárt alakú referenciaérték körlap esetén;
- `Hiba`: a referenciaértéktől vett abszolút eltérés;
- `Idő`: a mért futási idő másodpercben.

A teljes síkbeli és Monte Carlo változatoknál a referenciaérték `pi`, ezért a hiba:

```text
|pi - eredmeny|
```

A körlapos változatnál:

```text
|pi * (1 - e^(-R^2)) - eredmeny|
```

## Teljesítménymérési javaslat

A szekvenciális és párhuzamos verziók összehasonlításához ugyanazokat a bemeneti
paramétereket kell használni.

Példa méréssorozat teljes síkbeli rácsos integrálra:

```bash
./gauss2d_integral_seq 1000 1000
OMP_NUM_THREADS=2 ./gauss2d_integral_par 1000 1000
OMP_NUM_THREADS=4 ./gauss2d_integral_par 1000 1000
OMP_NUM_THREADS=8 ./gauss2d_integral_par 1000 1000
```

Mérhető mutatók:

```text
gyorsulas = T_szekvencialis / T_parhuzamos
hatekonysag = gyorsulas / szalak_szama
```

A pontos összehasonlításhoz érdemes minden mérést többször lefuttatni, majd az
átlagos futási időt használni, mert az operációs rendszer terhelése és a CPU
állapotai befolyásolhatják az eredményt.

## Pontossági megfontolások

- Nagyobb `N`, `M`, `N_r` és `N_phi` értékek általában pontosabb determinisztikus
  integrálást adnak, de növelik a futási időt.
- Monte Carlo esetén nagyobb mintaszám csökkenti a statisztikus hibát, de a hiba
  futásonként ingadozhat.
- A transzformált módszer a végtelen tartományt véges intervallumra képezi le, de
  a `tan` és `cos` függvények miatt a végpontok közelében numerikusan érzékenyebb
  lehet. A program középpontos mintavételezést használ, ezért nem értékeli ki pontosan
  a szinguláris végpontokat.
- A programok `long double` típust használnak az összegzéshez és a fontosabb
  lebegőpontos számításokhoz, ami javíthatja a numerikus pontosságot a sima `double`
  típushoz képest.

## Fontos implementációs részletek

### Időmérés

Minden program a `clock_gettime(CLOCK_MONOTONIC, ...)` hívást használja. Ez monoton
órával mér, tehát nem zavarja meg a rendszeridő kézi vagy hálózati módosítása.

### Matematikai konstansok

A kód saját `long double` pontosságú pi konstansokat definiál:

```c
#define M_PIl   3.141592653589793238462643383279502884L
#define M_PI_2l 1.570796326794896619231321691639751442L
```

Ez hordozhatóbb, mint arra hagyatkozni, hogy a rendszer matematikai fejlécei
definiálják-e az `M_PI` makrót.

### Redukció

A párhuzamos összegzésnél a `sum` változó közös írása versenyhelyzetet okozna.
Ezt az OpenMP `reduction(+:sum)` záradék oldja meg: minden szál külön részösszeget
vezet, majd a párhuzamos szakasz végén ezek összeadódnak.

### Véletlenszám-generátor

A szekvenciális Monte Carlo program:

```c
srand(time(NULL));
rand();
```

megoldást használ. A párhuzamos változat szálankénti maggal indítja a generátort:

```c
unsigned int seed = time(NULL) + omp_get_thread_num();
rand_r(&seed);
```

Ez csökkenti a globális állapotból adódó problémákat, de tudományos igényű Monte
Carlo szimulációhoz érdemes lehet jobb minőségű, reprodukálható generátorra váltani.

## Ismert korlátok és fejlesztési lehetőségek

- A bemeneti paraméterek ellenőrzése minimális. Negatív vagy nulla felosztásszám
  esetén a programok hibás eredményt vagy futási hibát adhatnak.
- Nincs közös build rendszer, például `Makefile` vagy CMake konfiguráció.
- A kimeneti szövegek a forrásfájlokban jelenleg hibás karakterkódolásúnak tűnnek
  bizonyos magyar ékezetes karaktereknél.
- A párhuzamos Monte Carlo futás nem teljesen reprodukálható, mert az idő alapján
  inicializálja a véletlen magokat.
- A `gauss2d` segédfüggvény néhány fájlban megmaradt, de a gyorsított körlapos
  változatban az `e^(-r^2)` közvetlen számolása miatt nincs aktívan használva.

Lehetséges továbbfejlesztések:

- `Makefile` hozzáadása az összes program egyszerű fordításához;
- paraméterellenőrzés és egyértelmű hibaüzenetek bevezetése;
- mérési script készítése több szál és problémaméret automatikus teszteléséhez;
- reprodukálható Monte Carlo futás opcionálisan megadható seed paraméterrel;
- eredmények CSV fájlba írása teljesítmény-összehasonlításhoz;
- jobb minőségű véletlenszám-generátor használata a Monte Carlo változatokban.

## Javasolt beadandó-struktúra

Ha a projekthez külön írásos beadandó vagy prezentáció készül, az alábbi fejezetek
használhatók:

1. Feladat ismertetése
2. Matematikai háttér: 2D Gauss-integrál és analitikus eredmények
3. Numerikus módszerek: középpontos téglalap-szabály és Monte Carlo
4. Párhuzamosítás OpenMP-vel
5. Fordítás és futtatás
6. Mérési eredmények
7. Gyorsulás és hatékonyság elemzése
8. Pontosság és hiba vizsgálata
9. Következtetések és továbbfejlesztési lehetőségek

