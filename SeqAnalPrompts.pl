#! usr/bin/perl


use strict;

my @files;
my @file_names;
my @headers;
my @sequences;
my @id;
my @inserts;
my @length_inserts;
my @frame_shift;

my $seq;
my $frame;
my $length_insert;

my $file_count;
my $file;
my $line;
my $sequence;
my $file_name;

my $const_frame_leader_all;

my $const_frame_constant_gamma;
my $const_frame_constant_kappa;
my $const_frame_constant_lamda;

my $constant_vector_gamma;
my $constant_vector_kappa;
my $constant_vector_lamda;



# vector sequence used for assessment of reading-frame. common leader fragment of all vectors, includes AgeI site
$const_frame_leader_all = "ctagtagcaactgcaaccggt"; 

# Vector sequences used for assessment of reading-frame
$const_frame_constant_gamma = "..gtcgaccaagggcccatc"; # Cheavy (hIgG1)
$const_frame_constant_kappa = "acggtggctgcaccatctgtc"; # Ckappa (hIgk)
$const_frame_constant_lamda = "gaggagcttcaagcc"; # Clambda (hIgl3), could also be shifted 5' for the MscI lambda vector


# Vectors sequences used for identification of cloning vectors. Located a couple of codons downstream of the 3' restriction site.
$constant_vector_gamma = "atcggtcttccccctggcaccctc";
$constant_vector_kappa = "ccatctgtcttcatcttcccgcca";
$constant_vector_lamda = "gtcactctgttccc";

$file_names[0]= "File";
$headers[0] = "Header";
$id[0] = "ID";
$length_inserts[0] = "Insert Length";
$inserts[0] = "Insert Sequence";
$frame_shift[0] = "Frame Shift";
$sequences[0] = "Sequence";
my $user_input;
print "enter directory location\n";
$user_input = <STDIN>;
chomp($user_input);
$file_count= 0;
@files = <$user_input/*>;
foreach $file (@files) {
	$file_name = $file;
	push (@file_names, $file_name);
	open(FILE, $file) or die "can't open file";
	while($line=<FILE>) {
		if ($line =~ /^>/) {
			$file_count++;
			chomp($line);
			push (@headers, $line);
			
			if(defined($sequence)) {
			
				$sequence = "";
			}		
		}
		else {
			chomp($line);
			$seq .= $line;
		}
	}
	$sequence= lc($seq);
	push (@sequences, $sequence);
	my ($seq, $conserved_region_sequence) = insert_id($sequence, $const_frame_leader_all, $constant_vector_gamma, $constant_vector_kappa, $constant_vector_lamda, $const_frame_constant_gamma, $const_frame_constant_kappa, $const_frame_constant_lamda);	

	close (FILE);
}




#print array of data
my $x = 0;
for ($x= 0; $x<$file_count+1; $x++){
	print " $file_names[$x] \t $headers[$x] \t $id[$x] \t $frame_shift[$x] \n"
}

my $x = 0;
for ($x= 0; $x<$file_count+1; $x++){
open FILE, ">>output.csv";
print FILE "$file_names[$x], $headers[$x], $id[$x], $length_inserts[$x], $inserts[$x], $frame_shift[$x]\n ";
}





sub insert_id {
	my($seq, $global, $id_gamma, $id_kappa, $id_lamda, $conserved_gamma, $conserved_kappa, $conserved_lamda)= @_;
	my $gamma = "gamma";
	my $kappa = "kappa";
	my $lamda = "lamda";
	my $fail = "not a sequence";
	my $start = "aa";
	my $end;
	if ($seq =~ /$global/){
		if ($seq =~ /$id_gamma/){
			push (@id, $gamma);
			$end = $conserved_gamma;
			$seq =~ /$start(.+)$end/;
			push (@inserts, $1);
				$length_insert = length ($1);
				push (@length_inserts, $length_insert);
				my $x= $length_insert / 3 ;
					if ($x =~ /\d\.(\d*)/) {
						if ($1 =~ /3+/){
							$frame= 1;
						}
						elsif ($1 =~ /6+/){
							$frame= 2;
						}
					}
					else{
						$frame = 0;
					}
			push (@frame_shift, $frame);
		}
		elsif ($seq =~ /$id_kappa/){
			push (@id, $kappa);
			#print "kappa\n";
			$end = $conserved_kappa;
			$seq =~ /$start(.+)$end/;
			push (@inserts, $1);
				$length_insert = length ($1);
				push (@length_inserts, $length_insert);
				my $x= $length_insert / 3 ;
					if ($x =~ /\d\.(\d*)/) {
						if ($1 =~ /3+/){
							$frame= 1;
						}
						elsif ($1 =~ /6+/){
							$frame= 2;
						}
					}
					else{
						$frame = 0;
					}
			push (@frame_shift, $frame);
		}
		else {
			push (@id, $lamda);
			$end = $conserved_lamda;
			$seq =~ /$start(.+)$end/;
			push (@inserts, $1);
				$length_insert = length ($1);
				push (@length_inserts, $length_insert);
				my $x= $length_insert / 3 ;
					if ($x =~ /\d\.(\d*)/) {
						if ($1 =~ /3+/){
							$frame= 1;
						}
						elsif ($1 =~ /6+/){
							$frame= 2;
						}
					}
					else{
						$frame = 0;
					}
			push (@frame_shift, $frame);
		}
	}
	else {
		push(@id, $fail);
		push(@inserts, $fail);
		push(@length_inserts, $fail);
		push(@frame_shift, $fail);
	}
}












#sub revcompDNA {                    				
# these are my subroutines that do the reverse complement. this one is for dna
#	my ($DNA) = @_;
#	my $rev = reverse ($DNA); 
#	$rev =~ tr/ACGTacgt/TGCAtgca/; 
#	return $rev;
#	}