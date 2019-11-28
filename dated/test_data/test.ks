CODONML (in paml version 4.8a, August 2014)  test_data/test.pal2nal
Model: several dN/dS ratios for branches for branches, 
Codon frequency model: F3x4
ns =   2  ls = 1377

Codon usage in sequences
----------------------------------------------------------------------
Phe TTT  31  31 | Ser TCT  38  38 | Tyr TAT  11  11 | Cys TGT   9   9
    TTC  16  16 |     TCC  14  14 |     TAC   7   7 |     TGC   6   6
Leu TTA  26  27 |     TCA  56  56 | *** TAA   0   0 | *** TGA   0   0
    TTG  47  47 |     TCG   6   6 |     TAG   0   0 | Trp TGG   8   8
----------------------------------------------------------------------
Leu CTT  39  39 | Pro CCT  38  38 | His CAT  18  18 | Arg CGT  16  16
    CTC  17  17 |     CCC   9   9 |     CAC   9   8 |     CGC   8   9
    CTA  13  12 |     CCA  25  25 | Gln CAA  46  46 |     CGA   4   4
    CTG  18  18 |     CCG  10  10 |     CAG  31  31 |     CGG   5   5
----------------------------------------------------------------------
Ile ATT  48  49 | Thr ACT  37  36 | Asn AAT  29  29 | Ser AGT  37  37
    ATC  13  13 |     ACC   6   6 |     AAC  11  11 |     AGC  11  11
    ATA  17  17 |     ACA  38  38 | Lys AAA  32  32 | Arg AGA  13  13
Met ATG  24  24 |     ACG   5   5 |     AAG  40  39 |     AGG   7   7
----------------------------------------------------------------------
Val GTT  30  30 | Ala GCT  45  45 | Asp GAT  62  62 | Gly GGT  23  23
    GTC  14  14 |     GCC   8   8 |     GAC  17  17 |     GGC  12  12
    GTA  21  21 |     GCA  36  36 | Glu GAA  53  54 |     GGA  26  26
    GTG  30  30 |     GCG   5   5 |     GAG  38  38 |     GGG   8   8
----------------------------------------------------------------------

Codon position x base (3x4) table for each sequence.

#1: Phacu.010G001200.1
position  1:    T:0.19971    C:0.22222    A:0.26725    G:0.31082
position  2:    T:0.29339    C:0.27306    A:0.29339    G:0.14016
position  3:    T:0.37110    C:0.12927    A:0.29484    G:0.20479
Average         T:0.28807    C:0.20818    A:0.28516    G:0.21859

#2: Phacu.WLD.010G000600.1
position  1:    T:0.20044    C:0.22150    A:0.26652    G:0.31155
position  2:    T:0.29412    C:0.27233    A:0.29267    G:0.14089
position  3:    T:0.37110    C:0.12927    A:0.29557    G:0.20407
Average         T:0.28855    C:0.20770    A:0.28492    G:0.21883

Sums of codon usage counts
------------------------------------------------------------------------------
Phe F TTT      62 | Ser S TCT      76 | Tyr Y TAT      22 | Cys C TGT      18
      TTC      32 |       TCC      28 |       TAC      14 |       TGC      12
Leu L TTA      53 |       TCA     112 | *** * TAA       0 | *** * TGA       0
      TTG      94 |       TCG      12 |       TAG       0 | Trp W TGG      16
------------------------------------------------------------------------------
Leu L CTT      78 | Pro P CCT      76 | His H CAT      36 | Arg R CGT      32
      CTC      34 |       CCC      18 |       CAC      17 |       CGC      17
      CTA      25 |       CCA      50 | Gln Q CAA      92 |       CGA       8
      CTG      36 |       CCG      20 |       CAG      62 |       CGG      10
------------------------------------------------------------------------------
Ile I ATT      97 | Thr T ACT      73 | Asn N AAT      58 | Ser S AGT      74
      ATC      26 |       ACC      12 |       AAC      22 |       AGC      22
      ATA      34 |       ACA      76 | Lys K AAA      64 | Arg R AGA      26
Met M ATG      48 |       ACG      10 |       AAG      79 |       AGG      14
------------------------------------------------------------------------------
Val V GTT      60 | Ala A GCT      90 | Asp D GAT     124 | Gly G GGT      46
      GTC      28 |       GCC      16 |       GAC      34 |       GGC      24
      GTA      42 |       GCA      72 | Glu E GAA     107 |       GGA      52
      GTG      60 |       GCG      10 |       GAG      76 |       GGG      16
------------------------------------------------------------------------------


Codon position x base (3x4) table, overall

position  1:    T:0.20007    C:0.22186    A:0.26688    G:0.31118
position  2:    T:0.29375    C:0.27269    A:0.29303    G:0.14052
position  3:    T:0.37110    C:0.12927    A:0.29521    G:0.20443
Average         T:0.28831    C:0.20794    A:0.28504    G:0.21871


Nei & Gojobori 1986. dN/dS (dN, dS)
(Note: This matrix is not used in later ML. analysis.
Use runmode = -2 for ML pairwise comparison.)

Phacu.010G001200.1  
Phacu.WLD.010G000600.1 0.4816 (0.0010 0.0020)

pairwise comparison, codon frequencies: F3x4.


2 (Phacu.WLD.010G000600.1) ... 1 (Phacu.010G001200.1)
lnL =-5546.878240
  0.00369 999.00000  0.64638

t= 0.0037  S=  1211.2  N=  2919.8  dN/dS=  0.6464  dN = 0.0011  dS = 0.0016
