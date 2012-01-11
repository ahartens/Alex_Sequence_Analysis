#! usr/bin/perl


use strict;
#use warnings;

my @files;
my @file_names;
my @headers;
my @sequences;
my @id;
my @inserts;
my @length_inserts;
my @frame_shift;
my @seq_lengths;
my @stops;
my $user_output_name;
my $rename;
my $file_type;

my $insert;
my $frame;
my $length_insert;
my $not_sequence;
my $fail;

my $file_count;
my $file;
my $line;
my $sequence;
my $file_name;
my $seq;
my $seq_length;
my $fail;
my $const_frame_leader_all;
my $const_frame_leader_all_revcomp;

my $const_frame_constant_gamma;
my $const_frame_constant_kappa;
my $const_frame_constant_lamda;

my $constant_vector_gamma;
my $constant_vector_kappa;
my $constant_vector_lamda;




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
$id[0] = "Constant";
$length_inserts[0] = "Insert Length";
$inserts[0] = "Insert Sequence";
$frame_shift[0] = "Frame";
$sequences[0] = "Sequence";
$stops[0] = "Stops";
$not_sequence= "not a sequence";
$fail = "";




my $directory_input;
my $quality_select;

print " \n \n \n \t:: Sequence Analysis ::\n \n  \t DIRECTIONS: This script should be placed in the directory containing the \n \t folders of sequences to be analyzed. Enter the name of the folder containing \n \t sequences to be analyzed when prompted. This will automatically be used \n \t as the name for the output file. If you wish to change the name of the\n \t output file, enter yes when prompted. The output file will be located in \n \t the original directory.";

print "\n \n \t ENTER THE FILE NAME OR FOLDER NAME CONTAINING SEQUENCES TO BE ANALYZED : \t ";
$directory_input = <STDIN>;
chomp($directory_input);

#print "\n \n \t WHICH QUALITY SEQUENCE FILE WOULD YOU LIKE TO USE?\n\n \t LOWEST (Q20) -- ENTER '1'\n\n \t MEDIUM Q30 -- ENTER '2' \n\n \t HIGH Q40 -- ENTER '3' \n \n \t "; 
#$quality_select = <STDIN>;
#chomp($quality_select);
#if ($quality_select == 1 ) {
#	$file_type = "seq";
#	}
#elsif ($quality_select == 2 ) {
#	$file_type = "seq.clipped";
#	}
#elsif ($quality_select == 3 ) {
#	$file_type = "q40.sequence.clipped";
#	}
#else {
#	print "ERROR: not valid input\n"
#	}

print "\n\n \t RENAME OUTPUT FILE? (yes/no) : ";
$rename= <STDIN>;


if ($rename =~ /yes/ ) {
	print "\n \n \t ENTER DESIRED NAME FOR THE OUTPUT FILE : \t ";
	$user_output_name = <STDIN>;
	chomp($user_output_name);
}

else {
	$user_output_name = $directory_input;
}



chomp($user_output_name);
chomp($directory_input);
$file_count= 0;
@files = <$directory_input/*>;

foreach $file (@files) {
#	if ($file =~ /\.$file_type$/){
	if ($file =~ /\.seq$/){

	$file_name = $file;
	push (@file_names, $file_name);
		$file_count++;

	open(FILE, $file) or die "can't open file";
	
	while($line=<FILE>) {

		
		if ($line =~ /^>/) {
			chomp($line);
			$line =~ s/\cM\cJ?//g;

			push (@headers, $line);
		
			if(defined($seq)) {
				$seq = "";
			}		
		}
		
		
		else{
			$line =~ s/\cM\cJ?//g;
			chomp($line);
			$seq .= $line;
		}
	}
	
	$sequence= lc($seq);
	

	
	if ($sequence =~ m/[^acgt]/){
		
		$sequence = $not_sequence;

		push(@sequences, $sequence);
		push(@headers, $sequence);
		push(@seq_lengths, $sequence);
		push(@id, $fail);
		push(@inserts, $fail);
		push(@length_inserts, $fail);
		push(@frame_shift, $fail);
		push(@stops, $fail);
	}
	
	
	else{
		$seq_length = length($sequence);
		push (@seq_lengths, $seq_length);
		push (@sequences, $sequence);
		insert_id($sequence, $const_frame_leader_all, $const_frame_leader_all_revcomp, $constant_vector_gamma, $constant_vector_kappa, $constant_vector_lamda, $const_frame_constant_gamma, $const_frame_constant_kappa, $const_frame_constant_lamda);	
	}
	
	
	close (FILE);
}



}
#print array of data
my $x = 0;
print "\n\n";
for ($x= 0; $x<$file_count+1; $x++){
	print " $headers[$x] \t $id[$x] \t $length_inserts[$x] \t $frame_shift[$x] \t $stops[$x] \n"
}

my $x = 0;
for ($x= 0; $x<$file_count+1; $x++){
open FILE, ">>$user_output_name.csv";
print FILE "$headers[$x], $id[$x], $length_inserts[$x], $frame_shift[$x], $stops[$x]\n ";
}















sub insert_id {
	
	my($seq, $global, $global_revcomp, $id_gamma, $id_kappa, $id_lamda, $conserved_gamma, $conserved_kappa, $conserved_lamda)= @_;
	my $gamma = "gamma";
	my $kappa = "kappa";
	my $lamda = "lamda";
	my $reverse = "reverse";
	my $start = $global;
	my $end;
	
	if ($seq =~ /$global/){
		
			
			if ($seq =~ /$id_gamma/) {
					
					push (@id, $gamma);
					$end = $conserved_gamma;
					$seq =~ /$start(.+)$end/;
					$insert = $1; 
					push (@inserts, $insert);
					$length_insert = length ($insert);
					push (@length_inserts, $length_insert);
					calculate_frame_shift($length_insert, $insert)
			
			} elsif ($seq =~ /$id_kappa/) {
					
					push (@id, $kappa);
					$end = $conserved_kappa;
					$seq =~ /$start(.+)$end/;
					$insert = $1;
					push (@inserts, $insert);
					$length_insert = length ($insert);
					push (@length_inserts, $length_insert);
					calculate_frame_shift($length_insert, $insert);
			
			} elsif ($seq =~ /$id_lamda/){
					
					push (@id, $lamda);
					$end = $conserved_lamda;
					$seq =~ /$start(.+)$end/;
					$insert = $1;
					push (@inserts, $insert);
					$length_insert = length ($insert);
					push (@length_inserts, $length_insert);
					calculate_frame_shift($length_insert, $insert);
			
			} else {
					
					push (@id, $fail);
					push(@inserts, $fail);
					push(@length_inserts, $fail);
					push(@frame_shift, $fail);
					push(@stops, $fail);
		
			}
			
	
	}
	
	
	
	
	
	elsif ($seq =~ /$global_revcomp/){
			
			revcompDNA($seq);
			push (@id, $reverse);
			print "reverse it\n";
	
	}
		
	
	
	
	
	
	else {
			
			push(@id, $fail);
			push(@inserts, $fail);
			push(@length_inserts, $fail);
			push(@frame_shift, $fail);
			push(@stops, $fail);
	
	}
	

}















sub calculate_frame_shift{
	#my $fail = "--";
	my $result;
	my ($value, $insert)= @_;
	my $result= $value % 3 ;
		if ($result  == 0) {	
			early_termination_check($insert);
		}
		else{
			push (@stops, $fail);
		}
		push (@frame_shift, $result);
}















sub early_termination_check {

	my ($string) = @_;

	# 'true' means early sequence terminates prematurely, 'false' means no stop codon found
	my $true = "TRUE";
	my $false = "FALSE";
	
	# define stop codons: 	
	my $char = 'tga';
	my $char_2 = 'taa';
	my $char_3 = 'tag'; 

	# offset set to zero for indexing function	
	my $offset = 0;
	my $offset_2 = 0;
	my $offset_3 = 0;
	
	my $check;
	my $check_2;
	my $check_3;

	
	
	
	# looks for index of stop codon. Does modulus function to see if it is in the reading frame. If it is an early termination sequence, returns 'true'
	my $result = index($string, $char, $offset);
	$check = $result % 3;
	
	
	if ($check == 0) {	
		
			push (@stops, $true);
	
	}
	
	
	
	
	else{
		
			while ($result != -1) {
				
					$check = $result % 3;
					
			
					if ($check == 0) {
					
							push (@stops, $true);
				
					}
   		
   		
   					$offset = $result + 1;
    				$result = index($string, $char, $offset);	
		
			}
  	
  	
  	
  	
  			if ($stops[$file_count] ne $true) {
  		
  					my $result_2 = index($string, $char_2, $offset_2);
  					$check_2 = $result_2 % 3;
  		
  		
  		
  					if ($check_2 == 0) {
  				
  							push (@stops, $true);
  				
  					}
  		
  		
  		
  					else{
  			
  							while ($result_2 != -1) {
  					
  									$check_2 = $result_2 % 3;
  				
  				
  									if ($check_2 == 0) {
  							
  											push (@stops, $true);
  						
  									}
  				
  				
  									$offset_2 = $result_2 +1;
  									$result_2 = index($string, $char_2, $offset_2);
  			
  							}
  		
  					}
  		
  		
  		
  		
  			if ($stops[$file_count] ne $true){
  			
  					my $result_3 = index($string, $char_3, $offset_3);
  					$check_3 = $result % 3;
  			
  			
  			
  					if ($check_3 == 0) {
  							
  							push (@stops, $true);
  			
  					}
  			
  			
  			
  			
  					else{
  				
  							while ($result_3 != -1) {
  							
  									$check_3 = $result_3 % 3;
  					
  					
  									if ($check_3 == 0) {
  											
  											push (@stops, $true);
  								
  									}
  				
  				
  									$offset_3 = $result_3 +1;
  									$result_3 = index($string, $char_3, $offset_3);
  				
  							}
  			
  					}
  		
  			}
  	
  	}
  	
  	
  	
  	
			if ($stops[$file_count] ne $true){
  		
  			push (@stops, $false);
  	
  	}


	}

}















sub revcompDNA {                    				
	my ($DNA) = @_;
	my $rev = reverse ($DNA); 
	$rev =~ tr/ACGTacgt/TGCAtgca/; 
	return $rev;
	}