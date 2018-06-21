# seqlen
Script to print a comparison of protein and transcript length using two files from infoseq.

## Usage

```
Usage: seqlen.R [options]


Options:
	-h, --help
		Show this help message and exit

	-t TRANSCRIPT, --transcript=TRANSCRIPT
		File 1 with infoseq transcript sequence lengths.

	-p PROTEIN, --protein=PROTEIN
		File 2 with infoseq protein sequence lengths.

	-x XLAB, --xlab=XLAB
		X label for plot. [default 'Transcript Length (bp)']

	-y YLAB, --ylab=YLAB
		Y label for plot. [default 'Protein Length (codons)']

	-o OUTPUT, --output=OUTPUT
		Output file for plot.
```

### Example plots

<img src="https://raw.githubusercontent.com/scastlara/seqlen/master/example_plots/lencomparison.png" width=500/>

<img src="https://raw.githubusercontent.com/scastlara/seqlen/master/example_plots/bigger.png" width=500/>
