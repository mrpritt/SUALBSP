# SUALBSP

Supplements to [New solution approaches for balancing assembly lines with setup times](https://doi.org/10.1016/j.cor.2025.107202) which solves the setup assembly line balancing and scheduling problem (SUALBSP).

## Code

You can find the code in the [src](src) directory. To compile use
```bash
cmake -B build -G Ninja
cmake --build build --config Release
```
and to run
```bash
./build/hs --iter 100 --riter 100 --dwlb -p ap --mprule --heavy --alg S_hoffmann --cbfs --uselm4 --maxtimemodel 12 --exactt 120 --maxnodescbfs 1500 --secondpass instance.alb
./build/hs --iter 100 --riter 100 --dwlb -p ap --mprule --heavy --alg S_hoffmann --cbfs --uselm4 --maxtimemodel 24 --exactt 240 --maxnodescbfs 4500 --secondpass instance.alb
```
(the first example runs CBFS in the low configuration, the second in the high configuration of the paper). Instances are in [ALB format](https://assembly-line-balancing.de/sualbsp/data-set-of-scholl-et-al-2013). You can call `hs` with no arguments to see commandline options.

## Experimental data

To reproduce the tables in the paper, you can run the GNU R script [tables.R](results/tables.R).

## How to cite

```bibtex
@Article{Pereira.Ritt/2025,
  author =   {Jordi Pereira and Marcus Ritt},
  title =    {New solution approaches for balancing assembly lines with setup times},
  journal =  {Computers and Operations Research},
  volume =   {183},
  number =   {107202},
  month =    nov, 
  doi =      {10.1016/j.cor.2025.107202},
  url =      {https://doi.org/10.1016/j.cor.2025.107202}
}
```
