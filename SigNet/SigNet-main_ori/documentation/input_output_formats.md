# Signet Input/Output Formatting

## Mutations Input

There are different input formats accepted by SigNet. The default one is the mutation count matrix:

### Counts

Each row should be a different sample and each column should be a trinucleotide mutation type.
The different mutation types should follows the conventional 96 mutation type classification that is based on the six substitution subtypes: C>A, C>G, C>T, T>A, T>C, and T>G, as well as the nucleotides immediately 5’ and 3’ to the mutation [[1]](#1). Therefore, the shape of the data should be nx96, where n is the number of samples. It should also include a header with the trinucleotide mutation type for each column (format: A[C>A]A), and the sample ID as the index for each row. 

| ID | A[C>A]A | A[C>A]C | A[C>A]G | A[C>A]T | C[C>A]A | C[C>A]C | C[C>A]G | C[C>A]T | G[C>A]A | G[C>A]C | G[C>A]G | G[C>A]T | T[C>A]A | T[C>A]C | T[C>A]G | T[C>A]T | A[C>G]A | A[C>G]C | A[C>G]G | A[C>G]T | C[C>G]A | C[C>G]C | C[C>G]G | C[C>G]T | G[C>G]A | G[C>G]C | G[C>G]G | G[C>G]T | T[C>G]A | T[C>G]C | T[C>G]G | T[C>G]T | A[C>T]A | A[C>T]C | A[C>T]G | A[C>T]T | C[C>T]A | C[C>T]C | C[C>T]G | C[C>T]T | G[C>T]A | G[C>T]C | G[C>T]G | G[C>T]T | T[C>T]A | T[C>T]C | T[C>T]G | T[C>T]T | A[T>A]A | A[T>A]C | A[T>A]G | A[T>A]T | C[T>A]A | C[T>A]C | C[T>A]G | C[T>A]T | G[T>A]A | G[T>A]C | G[T>A]G | G[T>A]T | T[T>A]A | T[T>A]C | T[T>A]G | T[T>A]T | A[T>C]A | A[T>C]C | A[T>C]G | A[T>C]T | C[T>C]A | C[T>C]C | C[T>C]G | C[T>C]T | G[T>C]A | G[T>C]C | G[T>C]G | G[T>C]T | T[T>C]A | T[T>C]C | T[T>C]G | T[T>C]T | A[T>G]A | A[T>G]C | A[T>G]G | A[T>G]T | C[T>G]A | C[T>G]C | C[T>G]G | C[T>G]T | G[T>G]A | G[T>G]C | G[T>G]G | G[T>G]T | T[T>G]A | T[T>G]C | T[T>G]G | T[T>G]T |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| sample1 | 41 | 23 | 5 | 16 | 27 | 31 | 9 | 25 | 17 | 15 | 7 | 19 | 15 | 15 | 3 | 38 | 85 | 7 | 12 | 19 | 14 | 5 | 5 | 10 | 16 | 12 | 6 | 15 | 13 | 5 | 5 | 89 | 123 | 21 | 38 | 33 | 41 | 40 | 49 | 52 | 43 | 41 | 50 | 43 | 32 | 44 | 34 | 81 | 26 | 11 | 15 | 21 | 9 | 63 | 14 | 22 | 13 | 14 | 47 | 22 | 22 | 11 | 13 | 20 | 67 | 28 | 29 | 46 | 29 | 39 | 21 | 27 | 20 | 21 | 75 | 27 | 25 | 36 | 16 | 71 | 14 | 13 | 10 | 19 | 6 | 18 | 8 | 10 | 8 | 7 | 15 | 7 | 4 | 12 | 16 | 26 |
| sample2 | 37 | 14 | 5 | 18 | 33 | 16 | 14 | 28 | 13 | 10 | 9 | 12 | 20 | 25 | 4 | 35 | 58 | 15 | 5 | 17 | 11 | 4 | 6 | 12 | 7 | 7 | 8 | 14 | 13 | 9 | 4 | 85 | 119 | 31 | 24 | 34 | 38 | 41 | 32 | 40 | 30 | 29 | 41 | 31 | 35 | 37 | 31 | 55 | 26 | 9 | 17 | 27 | 8 | 38 | 17 | 14 | 6 | 17 | 51 | 22 | 15 | 15 | 7 | 24 | 41 | 30 | 25 | 41 | 22 | 27 | 10 | 32 | 20 | 19 | 54 | 27 | 23 | 36 | 21 | 40 | 17 | 6 | 15 | 9 | 8 | 6 | 10 | 3 | 8 | 8 | 11 | 12 | 12 | 9 | 11 | 34 |
| sample3 | 48 | 24 | 5 | 21 | 41 | 28 | 12 | 35 | 18 | 11 | 9 | 15 | 25 | 17 | 11 | 32 | 88 | 13 | 2 | 21 | 10 | 3 | 1 | 10 | 7 | 7 | 6 | 10 | 20 | 12 | 2 | 78 | 150 | 37 | 37 | 58 | 49 | 47 | 47 | 37 | 38 | 40 | 57 | 55 | 39 | 42 | 31 | 76 | 20 | 20 | 15 | 20 | 15 | 57 | 15 | 16 | 7 | 13 | 55 | 17 | 17 | 17 | 8 | 31 | 67 | 18 | 18 | 39 | 22 | 50 | 18 | 33 | 24 | 12 | 77 | 36 | 30 | 36 | 14 | 68 | 17 | 7 | 19 | 13 | 6 | 13 | 8 | 12 | 2 | 9 | 19 | 12 | 11 | 14 | 17 | 54 |
| sample4 | 56 | 28 | 6 | 19 | 43 | 29 | 12 | 27 | 15 | 15 | 4 | 21 | 32 | 16 | 5 | 45 | 70 | 17 | 2 | 18 | 11 | 7 | 3 | 9 | 16 | 11 | 3 | 18 | 14 | 14 | 0 | 86 | 127 | 44 | 34 | 39 | 41 | 53 | 57 | 36 | 51 | 46 | 56 | 52 | 42 | 62 | 36 | 73 | 21 | 12 | 20 | 21 | 10 | 64 | 21 | 22 | 16 | 7 | 48 | 16 | 30 | 12 | 14 | 26 | 67 | 31 | 33 | 45 | 37 | 37 | 32 | 33 | 17 | 20 | 55 | 26 | 27 | 38 | 12 | 65 | 17 | 12 | 12 | 11 | 6 | 8 | 9 | 5 | 7 | 11 | 27 | 13 | 8 | 12 | 20 | 31 |
| sample5 | 67 | 32 | 5 | 27 | 62 | 36 | 11 | 52 | 32 | 25 | 13 | 30 | 47 | 39 | 8 | 67 | 83 | 12 | 6 | 14 | 11 | 7 | 0 | 13 | 12 | 8 | 10 | 17 | 9 | 13 | 5 | 125 | 162 | 53 | 56 | 56 | 48 | 71 | 60 | 56 | 62 | 70 | 87 | 67 | 52 | 69 | 46 | 91 | 30 | 19 | 24 | 39 | 16 | 55 | 29 | 19 | 21 | 17 | 57 | 22 | 30 | 18 | 18 | 32 | 95 | 41 | 36 | 59 | 36 | 58 | 30 | 44 | 33 | 27 | 90 | 34 | 41 | 36 | 21 | 57 | 19 | 7 | 22 | 11 | 7 | 16 | 9 | 12 | 8 | 11 | 15 | 17 | 8 | 10 | 20 | 46 |

### Bed



However, there are other formats that are also accepted. Specifying --input_format bed, the user can provide a bed file that, in each row, contains a mutation of the sample  or the samples. In this case, the columns expected in this file are: sample, chr, start, ref, alt.

It is important to specify the reference genome, otherwise the trinucleotide context mapping cannot be done.

### VCF

Finally, specifying --input_format vcf, the user can provide a VCF file or a list of VCF files with the mutations of each sample. In this case it is also important to specify the reference genome.

## Normalization Input



File containing the abundances in two columns. The first one should be the trinucleotide, and the second one should be the abundance of that trinucleotide. The default option for normalization is None which means that no specific normalization will be applied and the data will just be normally normalized by the total number of mutations in each sample.

An example of the abundances file look like this:

<div style="height: 50px; overflow-y: auto !important; border: 1px solid #ccc;">

<table>
    <thead>
        <tr>
            <th></th>
            <th></th>
        </tr>
    </thead>
    <tbody>
        <tr><td>AAA</td><td>698329</td></tr>
        <tr><td>AAC</td><td>436449</td></tr>
        <tr><td>AAG</td><td>737530</td></tr>
        <tr><td>AAT</td><td>416568</td></tr>
        <tr><td>CAA</td><td>600982</td></tr>
        <tr><td>CAC</td><td>519841</td></tr>
        <tr><td>CAG</td><td>909499</td></tr>
        <tr><td>CAT</td><td>490253</td></tr>
        <tr><td>GAA</td><td>747391</td></tr>
        <tr><td>GAC</td><td>481873</td></tr>
        <tr><td>GAG</td><td>738007</td></tr>
        <tr><td>GAT</td><td>453174</td></tr>
        <tr><td>TAA</td><td>243098</td></tr>
        <tr><td>TAC</td><td>307855</td></tr>
        <tr><td>TAG</td><td>174900</td></tr>
        <tr><td>TAT</td><td>285445</td></tr>
        <tr><td>ACA</td><td>570372</td></tr>
        <tr><td>ACC</td><td>530099</td></tr>
        <tr><td>ACG</td><td>216770</td></tr>
        <tr><td>ACT</td><td>425910</td></tr>
        <tr><td>CCA</td><td>761018</td></tr>
        <tr><td>CCC</td><td>664368</td></tr>
        <tr><td>CCG</td><td>338357</td></tr>
        <tr><td>CCT</td><td>686587</td></tr>
        <tr><td>GCA</td><td>612890</td></tr>
        <tr><td>GCC</td><td>698452</td></tr>
        <tr><td>GCG</td><td>289058</td></tr>
        <tr><td>GCT</td><td>617328</td></tr>
        <tr><td>TCA</td><td>567719</td></tr>
        <tr><td>TCC</td><td>557411</td></tr>
        <tr><td>TCG</td><td>190380</td></tr>
        <tr><td>TCT</td><td>495949</td></tr>
        <tr><td>AGA</td><td>797782</td></tr>
        <tr><td>AGC</td><td>680942</td></tr>
        <tr><td>AGG</td><td>648844</td></tr>
        <tr><td>AGT</td><td>427683</td></tr>
        <tr><td>CGA</td><td>231066</td></tr>
        <tr><td>CGC</td><td>289095</td></tr>
        <tr><td>CGG</td><td>336285</td></tr>
        <tr><td>CGT</td><td>178119</td></tr>
        <tr><td>GGA</td><td>759821</td></tr>
        <tr><td>GGC</td><td>640772</td></tr>
        <tr><td>GGG</td><td>540477</td></tr>
        <tr><td>GGT</td><td>359468</td></tr>
        <tr><td>TGA</td><td>627645</td></tr>
        <tr><td>TGC</td><td>608043</td></tr>
        <tr><td>TGG</td><td>776144</td></tr>
        <tr><td>TGT</td><td>469233</td></tr>
        <tr><td>ATA</td><td>258515</td></tr>
        <tr><td>ATC</td><td>419516</td></tr>
        <tr><td>ATG</td><td>586183</td></tr>
        <tr><td>ATT</td><td>378712</td></tr>
        <tr><td>CTA</td><td>287063</td></tr>
        <tr><td>CTC</td><td>557813</td></tr>
        <tr><td>CTG</td><td>895650</td></tr>
        <tr><td>CTT</td><td>485231</td></tr>
        <tr><td>GTA</td><td>208254</td></tr>
        <tr><td>GTC</td><td>346293</td></tr>
        <tr><td>GTG</td><td>572729</td></tr>
        <tr><td>GTT</td><td>306524</td></tr>
        <tr><td>TTA</td><td>256116</td></tr>
        <tr><td>TTC</td><td>487837</td></tr>
        <tr><td>TTG</td><td>426503</td></tr>
        <tr><td>TTT</td><td>438956</td></tr>
    </tbody>
</table>

</div>


### Signet Refitter Output

The output are 4 text files containing, for each sample:
1. The signature weight guesses.
2. Lower bound error bars.
3. Upper bound error bars.
4. The detector score.

Each row corresponds to a sample and each column corresponds to a different signature.

The weight's output would look like this: 

| | SBS1| SBS2| SBS3| SBS4| SBS5| SBS6| SBS7a| SBS7b| SBS7c| SBS7d| SBS8| SBS9| SBS10a| SBS10b| SBS11| SBS12| SBS13| SBS14| SBS15| SBS16| SBS17a| SBS17b| SBS18| SBS19| SBS20| SBS21| SBS22| SBS23| SBS24| SBS25| SBS26| SBS27| SBS28| SBS29| SBS30| SBS31| SBS32| SBS33| SBS34| SBS35| SBS36| SBS37| SBS38| SBS39| SBS40| SBS41| SBS42| SBS43| SBS44| SBS45| SBS46| SBS47| SBS48| SBS49| SBS50| SBS51| SBS52| SBS53| SBS54| SBS55| SBS56| SBS57| SBS58| SBS59| SBS60| SBS84| SBS85| SBS86| SBS87| SBS88| SBS89| SBS90| Unknown |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|sample1| 0.0| 0.011373657| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.09554387| 0.0| 0.028311478| 0.0| 0.025709331| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.026231365| 0.0| 0.102276586| 0.0| 0.0| 0.05594749| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.022387402| 0.0| 0.0| 0.05920576| 0.0| 0.0| 0.0| 0.20035994| 0.0| 0.0| 0.0| 0.0| 0.13489234| 0.0| 0.0| 0.0| 0.043502145| 0.0| 0.0| 0.0| 0.0| 0.0| 0.04187152| 0.0| 0.0| 0.0| 0.0| 0.12799643| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.024390638 |
| sample2| 0.0| 0.04407453| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.09683562| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.017493004| 0.010427136| 0.0| 0.12316588| 0.033200223| 0.0| 0.055034526| 0.0| 0.20925371| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.018737758| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.061912738| 0.03293125| 0.013076979| 0.0| 0.0| 0.068620875| 0.05564928| 0.0| 0.03608794| 0.0| 0.010005509| 0.0| 0.0| 0.09760389| 0.0| 0.0| 0.0| 0.0| 0.015889168 |
| sample3| 0.020623289| 0.0291368| 0.79501826| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.05595904| 0.0| 0.0| 0.0| 0.0| 0.0| 0.032755993| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.04497706| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.021529555 |
| sample4| 0.30297777| 0.013867504| 0.0| 0.0| 0.4331048| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.011840954| 0.0| 0.0| 0.0| 0.0| 0.0| 0.09353717| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.14116676| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0035050511 |
| sample5| 0.32464442| 0.0| 0.0| 0.0| 0.5292391| 0.04437277| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.01917128| 0.042479694| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.0| 0.040092707 |
