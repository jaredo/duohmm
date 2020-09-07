# duohmm

Documentation: 

[https://mathgen.stats.ox.ac.uk/genetics_software/duohmm/duohmm.html](https://mathgen.stats.ox.ac.uk/genetics_software/duohmm/duohmm.html)

Publication:

[A General Approach for Haplotype Phasing across the Full Spectrum of Relatedness](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004234)

Build:

```
export BOOST_ROOT=/my/boost/installation/ ##optional

git clone https://github.com/jaredo/duohmm.git
cd duohmm/
make -j 4
```

Docker:

```
make clean
docker build -t duohmm .
```
