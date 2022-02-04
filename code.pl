use strict;
use warnings;


#1
print "\n1:\n";
my $DNA1 = "ACGTTA";
print "DNA1 is $DNA1\n";
my $DNA2 = "CTAA";
print "DNA2 is $DNA2\n";
my $DNA3 = $DNA1.$DNA2;
print "Concatenated DNA3 is $DNA3\n";


#2
print "\n2:\n";
my $RNA = "ACUUUG";
print "Given RNA sequence is $RNA\n";
my $DNA = $RNA;
$DNA =~ s/U/T/g;
print "Corresponding DNA sequence is $DNA\n";


#3
print "\n3:\n";
my $DNA_3 = "GCTTAC";
print "DNA sequence is $DNA_3\n";
my $revcom_3 = reverse $DNA_3;
$revcom_3 =~ tr/ACGTacgt/TGCAtgca/;
print "Reverse Complement is $revcom_3\n";


#4
print "\n4:\n";
open (SEQUENCE, "<", "E:/human zinc finger protein.txt") or die $!;
my $var = <SEQUENCE>; 
print $var;
close (SEQUENCE);


#5
print "\n\n5:\n";
my $DNA_5 = "CGattgaGCAt\n";
print "DNA sequence is $DNA_5";
my $DNA_5new1 = $DNA_5;
$DNA_5new1 =~ tr/ACGT/acgt/;
print "DNA sequence in lowercase is $DNA_5new1";
my $DNA_5new2 = $DNA_5;
$DNA_5new2 =~ tr/acgt/ACGT/;
print "DNA sequence in uppercase is $DNA_5new2\n";
 

#6
print "\n6:\n"; 
print "Enter DNA sequence:";
my $sequence = <STDIN>;
chomp $sequence;
print $sequence;
my $length = length($sequence);
print "Length of the DNA sequence is $length";


#7
print "\n7:\n"; 
@array = ('A','T','G','C');
print "Array = @array \n";
$pop = pop(@array);
print "\n The last element is : \n $pop \n";
$shift = shift(@array);
print "\n The first element is : \n $shift \n";
print "\n Array after transformation : \n @array \n";
push(@array,$pop);
print "\n Array after push : \n @array \n";
unshift(@array,$shift);
print "\n Array after unshift : \n @array \n\n";


#8
print "\n8:\n"; 
print "Enter the DNA sequence: \n";
my $dna = <>;
chomp $dna;
my $revcom = reverse($dna);
$revcom =~ tr/ATGCatgc/TACGtacg/; 
if ($revcom eq $dna ) 
{
print "Given sequence is palindromic! \n";
print "Original sequence: $dna \n";
print "Reverse complement: $revcom \n";
}
else
{
print "Given sequence is not palindromic! \n";
}


#9
print "\n9:\n"; 
print "Enter the DNA sequence: \n";
my $dna = <>;
chomp $dna;
$dna =~ tr/atgc/ATGC/;
$len = length($dna);
@DNA = split( '', $dna );
foreach $base (@DNA) {
if ( $base eq 'A' ) {
++$count_of_A;
} elsif ( $base eq 'C' ) {
++$count_of_C;
} elsif ( $base eq 'G' ) {
++$count_of_G;
} elsif ( $base eq 'T' ) {
++$count_of_T;
} else {
print "Error! Unrecognized base detected: $base\n";
++$errors;
}
}
$A = $count_of_A/$len * 100;
$T = $count_of_T/$len * 100;
$G = $count_of_G/$len * 100;
$C = $count_of_C/$len * 100;
$u = $errors/$len * 100;
print "A = $A% \n";
print "C = $T%\n";
print "G = $G%\n";
print "T = $C%\n";
print "Unrecognized = $u%\n";


#10
print "\n10:\n"; 
print "Please type the filename of the protein sequence 
data: ";
$proteinfilename = <STDIN>;
chomp $proteinfilename;
unless ( open(PROTEINFILE, $proteinfilename) ) 
{
print "Cannot open file \"$proteinfilename\"\n\n";
exit;
}
@protein = <PROTEINFILE>;
close PROTEINFILE;
$protein = join( '', @protein);
$protein =~ s/\s//g;
do {
print "Enter a motif to search for: ";
$motif = <STDIN>;
chomp $motif;
if ( $protein =~ /$motif/ ) {
print "Found the motif!\n\n";
} else {
print "Couldnt find the motif.\n\n";
}
} until ( $motif =~ /^\s*$/ );


#11
print "\n11:\n";
print "Enter the DNA sequence: \n";
my $dna = <>;
chomp $dna;
$dna =~ tr/atgc/ATGC/;
for ($i = 0; $i<length ($dna)-1; $i = $i+3)
{
 my $triplet = substr ($dna,$i,3);
 push(@codon, $triplet);
}
print "Codons : @codon\n";


#12
print "\n12:\n";
open(FH, "sry.fasta");
my @DNA = <FH>;
shift(@DNA);
my $dna = join('',@DNA);
$dna =~ s/\s//g;
$len = length($dna);
print "Length of the sequence : $len \n";
my $motif = "ATG.*TAA";
my $cds;
if ($dna =~ m/($motif)/g)
{
$cds = $1;
print "Given CDS = $cds \n";
$len = length($cds);
print "Length of the predicted CDS: $len \n";
}
else
{
print "No CDS detected \n";
}


#13
print "\n13:\n";
my $file = "E:/human zinc finger protein.txt";

unless (open(FILE, $file)) {
print "Could not open file $file!\n";
exit;
}

print "Reading the file and printing line by line:\n";
while(my $protein = <FILE>) {
print $protein;
}
close FILE;


#14
print "\n14:\n";
open("sequence", "<", "ZNF148.fasta") || die "Error reading file.\n";
my @prot = ();
while (my $item = <sequence>) {
chomp($item);
push(@prot, $item);
}
close("sequence");
shift @prot;
$prot = join( '', @prot);
$prot =~ s/\s//g;
$prot = join "", @prot;
if ( $prot =~ /N.*TT/ ) {
print "I found it!\n\n";
} else {
print "I couldn\'t find it.\n\n";
}


#15
print "\n\n15:\n";
print "Enter DNA sequence:";

sub polylength {
	my $input = <STDIN>;
	($input) = @_;
	my($size);
	my @codon = ();
	for (my $i=0; $i < (length($input)-2) ; $i +=3) {
        
        my $codon = substr($input,$i,3);
	push(@codon, $codon);
	}

	$size = scalar(@codon);
	return $size
}

print polylength(my$input);
print " = the length of the polypeptide chain";


#16
print "\n\n16:\n";
my $dna = "ATGGTTGCACCACAACCGGTTGTGCTCTGCAGTTGA";
my $aminoacid = " ";
my $codon;

for (my $i=0; $i < (length($dna)-2); $i +=3) {
$codon = substr($dna,$i,3);
$aminoacid .= codon2aa($codon);
}

print "The DNA $dna\nis translated to the amino acid $aminoacid\n";

my $rna = $dna;
$rna =~ s/T/U/g;
print "Corresponding mRNA sequence is $rna\n";
exit;


sub codon2aa {
my($codon) = @_;
$codon = uc $codon;

my(%genetic_code) = (
"TCA" => "S",
"TCC" => "S",
"TCG" => "S",
"TCT" => "S",
"TTT" => "F",
"TTA" => "L",
"TTG" => "L",
"TAC" => "Y",
"TAT" => "Y",
"TAA" => "_",
"TAG" => "_",
"TGC" => "C",
"TGT" => "C",
"TGA" => "_",
"TGG" => "W",
"CTA" => "L",
"CTC" => "L",
"CTG" => "L",
"CTT" => "L",
"CCA" => "P",
"CCC" => "P",
"CCG" => "P",
"CCT" => "P",
"CAC" => "H",
"CAT" => "H",
"CAA" => "Q",
"CAG" => "Q",
"CGA" => "R",
"CGC" => "R",
"CGG" => "R",
"CGT" => "R",
"ATA" => "I",
"ATC" => "I",
"ATT" => "I",
"ATG" => "M",
"ACA" => "T",
"ACC" => "T",
"ACG" => "T",
"ACT" => "T",
"AAC" => "N",
"AAT" => "N",
"AAA" => "K",
"AAG" => "K",
"AGC" => "S",
"AGT" => "S",
"AGA" => "R",
"AGG" => "R",
"GTA" => "V",
"GTC" => "V",
"GTG" => "V",
"GTT" => "V",
"GCA" => "A",
"GCC" => "A",
"GCG" => "A",
"GCT" => "A",
"GAC" => "D",
"GAT" => "D",
"GAA" => "E",
"GAG" => "E",
"GGA" => "G",
"GGC" => "G",
"GGG" => "G",
"GGT" => "G",
);

if (exists $genetic_code{$codon}) {
return $genetic_code{$codon};
} else {
print STDERR "Bad codon \"$codon\"!!\n";
exit;
}
}

#17
print "\n\n17:\n";
sub processFASTA {
my ($file_name) = @_;
open("sequence", "<", $file_name) || die "Error reading file.\n"; 
my @dna = ();
while (my $item = <sequence>) {
if (substr($item, 0, 1) ne ">" and substr($item, 0, 1) ne "#") 
{ 
push(@dna, $item);
}
}
close("sequence");
my $dna = join "", @dna; 
$dna =~ s/ //g;
$dna =~ s/\n//g;
return $dna;
}
print "The DNA sequence in FASTA file
is:\n".processFASTA("sample_dna.fasta");


#18
print "\n18:\n";
srand(time|$$);
my @nucs = ("A", "T", "G", "C");
while (1) {
srand(time|$$);
my $randNuc = randomnucleotide();
print "Your guess: ";
my $input = <STDIN>;
chomp($input);

if ($input eq $randNuc) {
print "You are correct! The nucleotide I chose was ",
$randNuc, "\n\n";
exit;
} else {
print "You are wrong, I selected $randNuc.\n\n"
}
}

sub randomelement {
my(@array) = @_;
return $array[rand @array];
}

sub randomnucleotide {
my(@nucleotides) = ('A', 'C', 'G', 'T');
return randomelement(@nucleotides);
}
  

#19
print "\n19:\n";
my $DNA = 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT';
my $i;
my $mutant;
srand(time|$$);

$mutant = mutate($DNA);
print "\nMutate DNA Program\n\n";
print "\nHere is the original DNA:\n\n";
print "$DNA\n";
print "\nHere is the mutant DNA:\n\n";
print "$mutant\n";

for ($i=0 ; $i < 10 ; ++$i) {
$mutant = mutate($mutant);
print "The mutant DNA is \n$mutant\n\n";
print "\n\nNow, translating mutated DNA to protein in all 6 frames.\n\n";

for (my $i = 1; $i<=3; $i++) {
print "\n -------Reading Frame $i--------\n\n";
my $protein = translate_frame($mutant, $i);
print $protein, "\n\n"
}

my $revcom = revcom($mutant);
for (my $i = 1; $i<=3; $i++) {
print "\n -------Reading Frame ", $i+3, "--------\n\n";
my $protein = translate_frame($revcom, $i);
print $protein, "\n\n"
}

}
exit;

sub revcom {
my($dna) = @_;
my $revcom = reverse $dna;
$revcom =~ tr/ACGTacgt/TGCAtgca/;
return $revcom;
}

sub translate_frame {
my($seq, $start, $end) = @_;
my $protein;
unless($end) {
$end = length($seq);
}
return dna2peptide ( substr ( $seq, $start - 1, $end -$start +
1) );
}

sub dna2peptide {
my($dna) = @_;
my $protein = '';
for(my $i=0; $i < (length($dna) - 2) ; $i += 3) {
$protein .= codon2aa( substr($dna,$i,3) );
}
return $protein;
}

sub codon2aa {
my($codon) = @_;
$codon = uc $codon;
my(%genetic_code) = (
'TCA' => 'S', # Serine
'TCC' => 'S', # Serine
'TCG' => 'S', # Serine
'TCT' => 'S', # Serine
'TTC' => 'F', # Phenylalanine
'TTT' => 'F', # Phenylalanine
'TTA' => 'L', # Leucine
'TTG' => 'L', # Leucine
'TAC' => 'Y', # Tyrosine
'TAT' => 'Y', # Tyrosine
'TAA' => '_', # Stop
'TAG' => '_', # Stop
'TGC' => 'C', # Cysteine
'TGT' => 'C', # Cysteine
'TGA' => '_', # Stop
'TGG' => 'W', # Tryptophan
'CTA' => 'L', # Leucine
'CTC' => 'L', # Leucine
'CTG' => 'L', # Leucine
'CTT' => 'L', # Leucine
'CCA' => 'P', # Proline
'CCC' => 'P', # Proline
'CCG' => 'P', #Proline
'CCT' => 'P', #Proline
'CAC' => 'H', #Histidine
'CAT' => 'H', #Histidine
'CAA' => 'Q', #Glutamine
'CAG' => 'Q', #Glutamine
'CGA' => 'R', #Arginine
'CGC' => 'R', #Arginine
'CGG' => 'R', #Arginine
'CGT' => 'R', #Arginine
'ATA' => 'I', #Isoleucine
'ATC' => 'I', #Isoleucine
'ATT' => 'I', #Isoleucine
'ATG' => 'M', #Methionine
'ACA' => 'T', #Threonine
'ACC' => 'T', #Threonine
'ACG' => 'T', #Threonine
'ACT' => 'T', #Threonine
'AAC' => 'N', #Asparagine
'AAT' => 'N', #Asparagine
'AAA' => 'K', #Lysine
'AAG' => 'K', #Lysine
'AGC' => 'S', #Serine
'AGT' => 'S', #Serine
'AGA' => 'R', #Arginine
'AGG' => 'R', #Arginine
'GTA' => 'V', #Valine
'GTC' => 'V', #Valine
'GTG' => 'V', #Valine
'GTT' => 'V', #Valine
'GCA' => 'A', #Alanine
'GCC' => 'A', #Alanine
'GCG' => 'A', #Alanine
'GCT' => 'A', #Alanine
'GAC' => 'D', #Aspartic Acid
'GAT' => 'D', #Aspartic Acid
'GAA' => 'E', #Glutamic Acid
'GAG' => 'E', #Glutamic Acid
'GGA' => 'G', #Glycine
'GGC' => 'G', #Glycine
'GGG' => 'G', # Glycine
'GGT' => 'G', # Glycine
);

if(exists $genetic_code{$codon}) {
return $genetic_code{$codon};
}
else{
print STDERR "Bad codon \"$codon\"!!\n";
exit;
}
}

sub mutate {
my($dna) = @_;
my(@nucleotides) = ('A', 'C', 'G', 'T');
my($position) = randomposition($dna);
my($newbase) = randomnucleotide(@nucleotides);
substr($dna,$position,1,$newbase);
return $dna;
}

sub randomelement {
my(@array) = @_;
return $array[rand @array];
}

sub randomnucleotide {
my(@nucleotides) = ('A', 'C', 'G', 'T');
return randomelement(@nucleotides);
}

sub randomposition {
my($string) = @_;
return int rand length $string;
}


