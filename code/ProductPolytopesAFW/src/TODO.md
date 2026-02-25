# TODOS

Questa lista raccoglie i prossimi task da implementare o investigare nel codice (soprattutto per generazione istanze e per rendere piu' “informative” le differenze FW vs AFW).

## [ ] Edge-Biased Sampling Per `box_uniform` (`box_boundary_power`)
Perche':
`ConvexHullLMO` scorre tutti i punti; i punti interni e/o ridondanti sono puro costo. In `box_uniform` oggi campioniamo nell'interno del box, quindi molti punti non sono vertici dell'hull e non aiutano l'LMO.

Idea:
Mantenere `vertex_sampling: "box_uniform"` ma aggiungere un parametro che concentra le coordinate vicino a LB/UB, aumentando la probabilita' che i punti siano “utili” (piu' vicini a facce/edge/corner del box) e riducendo il numero di punti necessari.

Implementazione (sketch):
1. Aggiungere a YAML un parametro `box_boundary_power: p` con vincolo `p >= 1.0` (default `1.0` = comportamento attuale).
2. In `generate_polytope` quando `vertex_sampling=="box_uniform"`:
   - campiona `u ~ Uniform(0,1)`;
   - con prob. 1/2 usa `t = u^p`, altrimenti `t = 1 - u^p`;
   - mappa `x = LB + (UB-LB) * t` per ogni coordinata.
3. (Opzionale) Aggiungere `box_boundary_eps` per tenere i punti strettamente dentro: clamp `t` in `[eps, 1-eps]`.
4. Aggiornare `generate_filename` per includere `p` (es. `box-unif-p2.0`).

File potenzialmente da toccare:
`examples/config.yml`, `src/config.jl`, `src/polytope_generation.jl`, `src/utils.jl`.

## [ ] Regolarizzare La “Forma” Dei Box (`box_width_min_factor`)
Perche':
Con `sphere` il raggio usa `min(h[d])`, quindi basta una coordinata “stretta” per collassare la scala della sfera. Con `ellipsoid`, gli assi ereditano direttamente l'anisotropia dei box e possono risultare molto schiacciati. Entrambi i fenomeni peggiorano il conditioning.

Idea:
Controllare la distribuzione delle larghezze per coordinata dei box, imponendo una larghezza minima per ogni dimensione, cosi' `min(h)` non diventa minuscolo e i rapporti `max/min` restano gestibili.

Implementazione (sketch):
1. Aggiungere a YAML `box_width_min_factor: wmin` con vincolo `0 < wmin <= 1`.
2. In `generate_nonintersecting_bounds`, quando si costruisce `upper_bounds[d]`:
   - imporre `(UB-LB) >= wmin * stepsize` (ad es. campionare width in `[wmin*stepsize, stepsize]`).
3. (Opzionale) Stampare anche `max/min` delle half-lengths nei diagnostics (gia' stampiamo `max/min`).
4. Aggiornare filename con `wmin` se attivo.

File potenzialmente da toccare:
`examples/config.yml`, `src/config.jl`, `src/polytope_generation.jl`, `src/utils.jl`.

## [ ] Separazione “Non Monotona” Dei Box (Harder Non-Intersecting Per `sphere/ellipsoid`)
Perche':
La generazione attuale separa i box in TUTTE le coordinate (P2 sta “sopra” P1 coordinata-per-coordinata). Se in piu' i polytopi sono piccoli rispetto alla separazione, il problema non-intersecting puo' diventare banale (convergenza in 1-2 iterazioni).

Idea:
Generare box disgiunti ma senza imporre separazione coordinata-per-coordinata (p.es. separazione lungo una direzione random o in un sottoinsieme di coordinate), mantenendo la garanzia di disgiunzione.

Implementazione (sketch):
1. Aggiungere `bounds_separation: "monotone"|"random"` in YAML.
2. Per `"random"`, campionare centri dei box e garantire disgiunzione via vincolo su distanza L_infty (o su una coordinata casuale) piu' `margin`.
3. Verificare con assert economici: se i box sono disgiunti allora anche i polytopi (che stanno dentro) sono disgiunti.
4. Aggiornare filename con tag del metodo (`sep=mono|rand`).

File potenzialmente da toccare:
`examples/config.yml`, `src/config.jl`, `src/polytope_generation.jl`, `src/utils.jl`.

## [ ] Preconditioning Globale (none/diag/whiten)
Perche':
Anisotropie e scale molto diverse tra coordinate peggiorano diametri/costanti geometriche; spesso AFW/FW si comportano meglio dopo centratura e scaling. Importante: la trasformazione deve essere la stessa per TUTTI i polytopi, altrimenti distruggi la proprieta' di intersezione/non-intersezione.

Idea:
Applicare una trasformazione affine globale ai vertici (e, se presenti, ai box bounds/LMOs strutturati):
1. `diag`: sottrai media globale e dividi per std (o range) per coordinata.
2. `whiten`: whitening (PCA) con regolarizzazione.

Implementazione (sketch):
1. Aggiungere a YAML: `preconditioning: "none"|"diag"|"whiten"`, e parametri tipo `whiten_eps`.
2. Implementare helper in `polytope_utils.jl` per calcolare `(mu, scale)` o matrice `W`.
3. In `generate_polytopes`, dopo aver generato e shiftato i polytopi, trasformare tutti i vertici.
4. Se/Quando introduciamo BoxLMO/HypercubeLMO, trasformare anche i bounds coerentemente (solo per `diag`; per `whiten` serve gestire un LMO non-box oppure evitare whiten quando ci sono box).
5. Aggiornare filename con tag (`pc=diag|whiten`).

File potenzialmente da toccare:
`examples/config.yml`, `src/config.jl`, `src/polytope_utils.jl`, `src/polytope_generation.jl`, `src/utils.jl`.

## [ ] Anchor `p1_random` “Piu' Profondo” (Dirichlet Alpha > 1)
Perche':
`p1_random` oggi e' una convex combination casuale; se i pesi sono troppo “sparsi”/vicini a 0, l'anchor puo' finire vicino a facce (e quindi reintrodurre problemi di conditioning locale).

Idea:
Sostituire i pesi random uniformi con un campionamento Dirichlet con parametro `alpha > 1` per spingere verso pesi piu' uniformi (anchor piu' interno).

Implementazione (sketch):
1. Aggiungere in YAML: `p1_random_dirichlet_alpha: 1.0` (default 1 = comportamento attuale).
2. Implementare Dirichlet senza dipendenze: campiona `w_i ~ Gamma(alpha,1)`, poi normalizza `w / sum(w)`.
3. Usare questi pesi in `random_convex_combination`.
4. Aggiornare filename con `diralpha`.

File potenzialmente da toccare:
`examples/config.yml`, `src/config.jl`, `src/polytope_utils.jl`, `src/utils.jl`.

## [ ] Istanza “Vertex–Facet” (Suggerimento Sebastian) Con LMOs Strutturati
Perche':
Costruire istanze che discriminano FW vs AFW: FW identifica una faccia grande solo asintoticamente, AFW puo' identificare in tempo finito grazie alle away steps (face/active-set identification).

Idea:
Per `k=2`:
1. P1 = box/cubo (LMO legacy, es. `BoxLMO` o `ZeroOneHypercubeLMO`).
2. P2 = simplex (convex hull di 3 punti) con un vertice che tocca il relint di un facet di P1 e gli altri due fuori dal box lungo la normale del facet, in modo che l'intersezione sia esattamente un punto.

Implementazione (sketch):
1. Estendere config per scegliere “tipo” di polytope per ogni blocco: `polytope_types: ["box","cvxhull"]` o simile.
2. Creare un generatore dedicato (nuova “family”) che costruisce direttamente:
   - bounds per P1 (BoxLMO),
   - vertici del simplex P2 (ConvexHullLMO) che realizzano il contatto vertex–facet.
3. Collegare questa family agli script di esempio (custom_full_warmup, slurm loop).
4. Aggiornare plotting/filenames per distinguere questa family dalle random-hull.

File potenzialmente da toccare:
`examples/config.yml`, `src/config.jl`, `src/lmo_utils.jl`, `src/polytope_generation.jl`, `src/utils.jl`, `examples/compute_intersection_custom_full_warmup.jl`.

## [ ] Supportare Polytopi “Legacy” E Misti (Box, Simplex, L1, ...)
Perche':
Permette di controllare meglio la geometria (width/diameter), evitare hull enormi, e costruire casi dove la teoria e' piu' interpretabile.

Idea:
Permettere che alcuni blocchi usino LMOs strutturati (BoxLMO, SimplexLMO, L1-ball, ecc.) e altri restino ConvexHullLMO.

Implementazione (sketch):
1. Estendere `Config` con una lista di tipi e parametri per polytope.
2. In `create_lmos`, creare il corrispondente LMO.
3. In `generate_polytopes`, se un blocco e' strutturato, non generare una nube di punti per quel blocco.

File potenzialmente da toccare:
`examples/config.yml`, `src/config.jl`, `src/lmo_utils.jl`, `src/polytope_generation.jl`, `src/ProductPolytopesAFW.jl`.

## [ ] “Pruning” Dei Vertici (Subsampling Estremo Per Accelerare ConvexHullLMO)
Perche':
Anche con edge-bias, `n_points` puo' restare enorme. Possiamo ridurre il set di atomi mantenendo (approssimativamente) l'hull, prendendo solo punti “estremi” lungo molte direzioni.

Idea:
Dato un cloud `V`, scegli M direzioni random `c_j` e conserva `argmin_v <c_j,v>` e/o `argmax_v <c_j,v>`; unisci e deduplica. Questo tende a selezionare vertici (o quasi-vertici) e riduce drasticamente cardinalita'.

Implementazione (sketch):
1. Aggiungere in YAML: `vertex_pruning: "none"|"random_directions"`, `num_directions: M`.
2. Implementare in `polytope_utils.jl` una funzione che restituisce un sottoinsieme di righe di `vertices`.
3. Chiamarla in `generate_polytope` (o subito dopo) prima di costruire l'LMO.
4. Aggiornare filename con `prune=M`.

File potenzialmente da toccare:
`examples/config.yml`, `src/config.jl`, `src/polytope_utils.jl`, `src/polytope_generation.jl`, `src/utils.jl`.

## [ ] Verifica Opzionale Intersezione (Cheap Check Via FW)
Perche':
Alcune configurazioni (es. `p1_p2_segment` con `anchor_t > 1`) possono rompere la garanzia di intersezione; inoltre errori numerici/bug sono sempre possibili.

Idea:
Opzionalmente, dopo lo shifting, lanciare un FW “opt” breve (o con `max_iterations_opt`) per stimare la distanza e assertare:
1. intersecting: `primal <= eps_intersection`,
2. non-intersecting: `primal > eps_intersection`.

Implementazione (sketch):
1. Aggiungere `verify_intersection: true|false` e `intersection_eps` a YAML.
2. In `generate_polytopes`, fare il check solo se richiesto (per non rallentare i batch).

File potenzialmente da toccare:
`examples/config.yml`, `src/config.jl`, `src/polytope_generation.jl`.

## [ ] Migliorare Filenames/Metadata Per Riproducibilita'
Perche':
Molti parametri influenzano fortemente i plot (es. `n_points` effettivo, sampling mode, preconditioning). Senza un nome/metadata completo e' facile confondersi e non riprodurre i risultati.

Idea:
1. Includere nel filename almeno un riassunto di `n_points` effettivo: `npmin/npmax/npmean` o la lista se corta.
2. Salvare nel `.jld2` anche la configurazione completa e le statistiche geometriche (shrink, max/min, ecc.).

Implementazione (sketch):
1. Aggiornare `generate_filename` e `save_polytopes`.
2. Aggiornare `plot_from_logs.jl` / `plotting_utils.jl` per leggere metadata quando utile.

File potenzialmente da toccare:
`src/utils.jl`, `src/polytope_generation.jl`, `src/plotting_utils.jl`, `examples/plot_from_logs.jl`.

## [ ] Rendere I Diagnostics (shrink stats) Opzionali E Consistenti
Perche':
Le stampe sono utili in locale ma rumorose su SLURM; inoltre vogliamo lo stesso formato per tutti i sampling mode (box_uniform / sphere / ellipsoid) e per future opzioni.

Idea:
Aggiungere un flag di logging e centralizzare il calcolo delle statistiche in una funzione utility.

Implementazione (sketch):
1. `print_geometry_stats: true|false` in YAML.
2. Spostare la logica di stampa in `polytope_utils.jl`.

File potenzialmente da toccare:
`examples/config.yml`, `src/config.jl`, `src/polytope_utils.jl`, `src/polytope_generation.jl`.

## [ ] Fix `log_times`: Tot Time Non Deve Sommare Tempi Cumulativi
Perche':
`log_times` fa `sum(iter -> iter[5], traj)`; se `iter[5]` e' “tempo cumulativo” (tipico), la somma esplode (vedi totali assurdi nei log).

Idea:
Usare l'ultimo timestamp come totale (se cumulativo) oppure differenziare (delta) tra iterazioni.

Implementazione (sketch):
1. Capire la semantica di `iter[5]` nelle trajectories (cumulativo vs delta).
2. Se cumulativo: `tot_time = traj[end][5]`.
3. Se delta: lasciare `sum`.
4. Aggiornare stampa e CSV.

File potenzialmente da toccare:
`src/utils.jl`, `examples/compute_intersection_custom_full_warmup.jl`.

## [ ] Gap Plottati: Primal Gap vs FW Gap (Globale vs Locale)
Perche':
In alcuni casi abbiamo visto FW gap < primal; possibile mismatch tra gap globale e gap di blocco (per block-coordinate) o mismatch nel plotting (uso/non uso di `f*`).

Idea:
1. Verificare cosa ritorna FrankWolfe.jl come “gap” per i metodi block-coordinate.
2. Se serve, calcolare un gap globale (max su blocchi) per plotting, almeno in modalita' debug.
3. Assicurarsi che nei plot si usi sempre `f(x_k) - f*` dove `f*` e' stimato con run “opt”.

Implementazione (sketch):
1. Ispezionare come vengono costruite le trajectories in `product_algorithms.jl`.
2. Eventualmente loggare sia gap locale sia gap globale.
3. Aggiornare `plotting_utils.jl`.

File potenzialmente da toccare:
`src/product_algorithms.jl`, `src/plotting_utils.jl`, `src/utils.jl`, `examples/compute_intersection_custom_full_warmup.jl`.
