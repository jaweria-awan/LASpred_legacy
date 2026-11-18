# LASpred (legacy, 2011)

The original Java implementation of LASpred, a lysine acetylation site scorer,
written in 2011 as a NetBeans project. It is kept here **as a reference, not
for active use**. The maintained version is the Rust port:
**https://github.com/jaweria-awan/LASpred**. The port uses the same model and
fixes the bugs listed below.

`src/finalcode/Main.java` is unchanged from 2011.

## What it does

The program reads a protein FASTA file and finds every lysine (K) with a
21-residue window around it. It scores each window as a weighted sum of nine
features: seven amino-acid property scales, a position frequency profile built
from `InputFileForHumanPosi.txt` (2,000 known human sites), and a dipeptide
table. Sites scoring above 212.368 are written to the output.

## Contents

| File | Purpose |
|---|---|
| `src/finalcode/Main.java` | The whole program |
| `InputFileForHumanPosi.txt` | Training windows: 2,000 human acetylation sites, 21 residues, K in the middle |
| `example.fasta` | Small test input: histone H3.1, histone H4, ubiquitin and one messy record |
| `run.sh` | Compiles and runs the code in Docker |
| `build.xml`, `manifest.mf`, `nbproject/` | NetBeans/Ant project files |

## How to run

The code has three Windows paths built in (`E:/e-coli15pos.txt` for input,
`E:/outpute-coli.txt` for output, `E:/InputFileForHumanPosi.txt` for training).
`run.sh` replaces them in a temporary copy, so `Main.java` stays untouched.

Requirements: Docker. No Java installation is needed.

```sh
git clone git@github.com:jaweria-awan/LASpred_legacy.git && cd LASpred_legacy
./run.sh example.fasta output.txt
```

Use your own proteins by putting the FASTA file in this folder:

```sh
./run.sh my_proteins.fasta my_output.txt
```

### Output

`output.txt` repeats each FASTA header and sequence, followed by one entry per
site above the threshold:

```
Position Of K14      KQTARKSTGGKAPRKQLATKA      283.4109196
```

- **Position:** 0-based (add 1 for the usual protein numbering).
- **Window:** the 21 residues around the K.
- **Number:** the integrated score.

### Without Docker

Edit the three `E:/...` paths in `Main.java`, then compile and run with any
JDK:

```sh
javac src/finalcode/Main.java -d out
java -cp out finalcode.Main
```

The NetBeans/Ant build (`ant jar`) is configured for Java 1.5
(`javac.source=1.5` in `nbproject/project.properties`). Modern JDKs no longer
support that setting, so use the `javac` commands above instead.

## Known problems

These are all fixed in the Rust version.

- **Threshold:** 212.368 is far below the scores of real windows (about 284 for
  known sites, about 276 for random ones), so almost every lysine is reported.
- **Off-by-one:** the first and last lysine with a full 21-residue window are
  never scored.
- **Input must start with `>`:** otherwise the program loops forever. `run.sh`
  checks for this first.
- **Lowercase and unknown residues:** lowercase sequences find no lysines at all
  (the `messy_example` record in `example.fasta` shows this). Residues such as
  X contribute 0 to the score instead of being flagged.
- **Missing dipeptides:** CW, MW and WM have no value in the dimer table and
  contribute 0.
- **Output format:** sites below the threshold print without a line break, so
  entries can run together on one line.
- **Circular evaluation:** the position profile is built from the training
  sites, so scoring those same sites overstates accuracy.
- **Compiled files:** the `.class`/`.jar` files from 2011 are not included, as
  they were built before the last source change.
