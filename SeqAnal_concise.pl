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
my @seq_lengths;



my $frame;
my $length_insert;
#my $fail;

my $file_count;
my $file;
my $line;
my $sequence;
my $file_name;
my $seq;
my $sequence;
my $seq_length;
my $fail;
my $skip;
my $const_frame_leader_all;
my $const_frame_leader_all_revcomp;
my $conserved_region_sequence;

my $const_frame_constant_gamma;
my $const_frame_constant_kappa;
my $const_frame_constant_lamda;

my $constant_vector_gamma;
my $constant_vector_kappa;
my $constant_vector_lamda;


#	my $fail = "--";


# Vector sequence used for assessment of reading-frame. common leader fragment of all vectors, includes AgeI site
$const_frame_leader_all = "ctagtagcaactgcaaccggt"; 
$const_frame_leader_all_revcomp = "accggttgcagttgctactag";

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




$file_count= 0;
@files = <analysis_example/*>;
foreach $file (@files) {
	$file_name = $file;
	push (@file_names, $file_name);
	open(FILE, $file) or die "can't open file";
	while($line=<FILE>) {
		if ($line =~ /^>/) {
			$file_count++;
			chomp($line);
			print "$line\n";
			push (@headers, $line);
			
			if(defined($seq)) {
			
				$seq = "";
			}		
		}
		else {
			chomp($line);
			$seq .= $line;
		}
	}
	$sequence= lc($seq);
	$seq_length = length($sequence);
	push (@seq_lengths, $seq_length);
	push (@sequences, $sequence);
	my ($global, $conserved_region_sequence) = insert_id($sequence, $const_frame_leader_all, $const_frame_leader_all_revcomp, $constant_vector_gamma, $constant_vector_kappa, $constant_vector_lamda, $const_frame_constant_gamma, $const_frame_constant_kappa, $const_frame_constant_lamda, $fail);	
	calculate_frame_shift($sequence, $global, $conserved_region_sequence, $fail);
	close (FILE);
}




#print array of data
my $x = 0;
for ($x= 0; $x<$file_count+1; $x++){
	print "$headers[$x] \t $id[$x] \t $seq_lengths[$x] \t $length_inserts[$x] \t $frame_shift[$x] \n"
}

my $x = 0;
for ($x= 0; $x<$file_count+1; $x++){
open FILE, ">>output.csv";
print FILE "$file_names[$x], $headers[$x], $id[$x], $length_inserts[$x], $frame_shift[$x]\n ";
}





sub insert_id {
	my($seq, $global, $global_revcomp, $id_gamma, $id_kappa, $id_lamda, $conserved_gamma, $conserved_kappa, $conserved_lamda, $fail)= @_;
	my $gamma = "gamma";
	my $kappa = "kappa";
	my $lamda = "lamda";
	my $start = $global;
	my $end;
	if ($seq =~ /$global/){
		if ($seq =~ /$id_gamma/){
			push (@id, $gamma);
			return ($global, $conserved_gamma);
		}
		elsif ($seq =~ /$id_kappa/){
			push (@id, $kappa);
			return ($global, $conserved_kappa);
		}
		else {
			push (@id, $lamda);
			return ($global, $conserved_lamda);
		}
	}
	elsif ($seq =~ /$global_revcomp/){
			return ($global_revcomp, $conserved_lamda);
		}
	else{
			$fail = "--";
			push(@id, $fail);
			return($fail);
	}
}





sub calculate_frame_shift{
	my ($seq, $start, $end) = @_;
	my $result;
	my $value;
	if ( defined($fail)){
		push(@inserts, $fail);
		push(@length_inserts, $fail);
		push(@frame_shift, $fail);
	}
	else {
		#print "$seq \n $start \t $end \n\n";
			
		$seq =~ /$start(.+)$end/;
			push (@inserts, $1);
		$length_insert = length ($1);
			push (@length_inserts, $length_insert);
		print "$length_insert\n";
		my $value= $length_insert / 3 ;
		if ($value =~ /\d\.(\d*)/) {
			if ($1 =~ /3+/){
				$result= 1;
				push (@frame_shift, $result);
			}
			elsif ($1 =~ /6+/){
				$result= 2;
				push (@frame_shift, $result);
			}
		}
		else{
			$result = 0;
			push (@frame_shift, $result);
		}
	
		}
	
		my $fail = "";

}








#sub revcompDNA {                    				
# these are my subroutines that do the reverse complement. this one is for dna
#	my ($DNA) = @_;
#	my $rev = reverse ($DNA); 
#	$rev =~ tr/ACGTacgt/TGCAtgca/; 
#	return $rev;
#	}