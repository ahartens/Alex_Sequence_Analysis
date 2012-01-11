#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

my @files;
my @headers;
my @sequences;
my @id;
my @inserts;
my @length_inserts;
my @frame_shift;
my @seq_lengths;
my $user_output_name;
my $rename;
my $file_type;
my $leader;
my %files;
my $file_name;
my $infos;
my @stop;
my $termination;
my $stop;
my $frame_shift;
my $insert;
my $frame;
my $length_insert;
my $not_sequence;
my $fail;
my $result;

my $id;
my $file_count;
my $file;
my $line;
my $sequence;
my $seq;
my $seq_length;

# Vector sequence used for assessment of reading-frame. common leader fragment of all vectors, includes AgeI site
my $const_frame_leader_all = "ctagtagcaactgcaaccggt"; 
my $const_frame_leader_all_revcomp = "accggttgcagttgctactag";

# Vector sequences used for assessment of reading-frame
my $const_frame_constant_gamma = "..gtcgaccaagggcccatc"; # Cheavy (hIgG1)
my $const_frame_constant_kappa = "acggtggctgcaccatctgtc"; # Ckappa (hIgk)
my $const_frame_constant_lamda = "gaggagcttcaagcc"; # Clambda (hIgl3), could also be shifted 5' for the MscI lambda vector

# Vectors sequences used for identification of cloning vectors. Located a couple of codons downstream of the 3' restriction site.
my $constant_vector_gamma = "atcggtcttccccctggcaccctc";
my $constant_vector_kappa = "ccatctgtcttcatcttcccgcca";
my $constant_vector_lamda = "gtcactctgttccc";


$not_sequence= "not a sequence";
$fail = "--";

my $directory_input;



if (defined($ARGV[0])) {
    $directory_input= $ARGV[0];
    chomp($directory_input);
	@files = <$directory_input/*>;
	
	$file_count= 0;

	foreach $file (@files) {
    	if ($file =~ /\.seq$/){
    	    $file_name = $file;
            $file_count++;
            $termination = "";

        	open(FILE, $file) or die "can't open file";
        	    while($line=<FILE>) {
            		if ($line =~ /^>/) {
                	    chomp($line);
                		$line =~ s/\cM\cJ?//g;
               			$files{$file_name}{header} = $line;
                    if(defined($seq)) {
                        $seq = "";
                    }
                } else{
                    $line =~ s/\cM\cJ?//g;
                    chomp($line);
                    $seq .= $line;
                }
            }
            $sequence= lc($seq);

            if ($sequence =~ m/[^acgt]/){

                $sequence = $not_sequence;
                $files{$file_name}{header} = $fail;                
                $files{$file_name}{sequence} = $fail;
                $files{$file_name}{seq_length} = $fail;
            	$files{$file_name}{leader} = $fail;
            	$files{$file_name}{constant} = $fail;
            	$files{$file_name}{insert} = $fail;
            	$files{$file_name}{length_insert} = $fail;
            	$files{$file_name}{frame} = $fail;
        	    $files{$file_name}{stop} = $fail;
        	    $files{$file_name}{file_count} = $file_count;


            } else {
                $seq_length = length($sequence);
                $files{$file_name}{sequence} = $sequence;
                $files{$file_name}{seq_length} = $seq_length;

                insert_id($sequence, $const_frame_leader_all, $const_frame_leader_all_revcomp, $constant_vector_gamma, $constant_vector_kappa, $constant_vector_lamda, $const_frame_constant_gamma, $const_frame_constant_kappa, $const_frame_constant_lamda);	
                
                $files{$file_name}{leader} = $leader;        
                $files{$file_name}{constant} = $id;
                $files{$file_name}{insert} = $insert;
                $files{$file_name}{length_insert} = $length_insert;
                $files{$file_name}{frame} = $result;
                $files{$file_name}{stop} = $termination;
                $files{$file_name}{file_count} = $file_count;

            }
            close (FILE);
        }
    }

    if ($file_count > 0){

        #  PRINT TERMINAL  #
        print " file count \t file name \t header \t sequence length \t constant \t insert length \t frame \t stop \n";
        for $file_name ( sort keys %files ) {
            print "$files{$file_name}{file_count} \t $file_name \t $files{$file_name}{header} \t $files{$file_name}{seq_length} \t  $files{$file_name}{constant} \t $files{$file_name}{length_insert} \t $files{$file_name}{frame} \t $files{$file_name}{stop} \n";
        }

        #  PRINT CSV  #
        open FILE, ">>$directory_input.csv";   
        print FILE " Filename, Constant, Leader, Length, Frame, Stop \n";

        for $file_name ( sort keys %files ) {
            print FILE "$file_name, $files{$file_name}{constant}, $files{$file_name}{leader}, $files{$file_name}{length_insert}, $files{$file_name}{frame}, $files{$file_name}{stop} \n";
        }
     } else{
	    print "no files\n";
	}



} else {
	print "\n \n \n ERROR: enter the name of a directory containing sequences to be analyzed\n\n\n";
}


# ========== SUBROUTINES ========== #

sub insert_id {
    my($seq, $global, $global_revcomp, $id_gamma, $id_kappa, $id_lamda, $conserved_gamma, $conserved_kappa, $conserved_lamda)= @_;
    my $gamma = "gamma";
    my $kappa = "kappa";
    my $lamda = "lamda";
    my $reverse = "reverse";
    my $start = $global;
    my $end;

    if ($seq =~ /$global/) {
    		$leader = "TRUE";
            if ($seq =~ /$id_gamma/) {
					$id = "gamma";
                    $end = $conserved_gamma;
                    $seq =~ /$start(.+)$end/;
                    $insert = $1;                     
                    $length_insert = length ($insert);
                    calculate_frame_shift($length_insert, $insert);
            } elsif ($seq =~ /$id_kappa/) {
					$id = "kappa";
                    $end = $conserved_kappa;
                    $seq =~ /$start(.+)$end/;
                    $insert = $1;
                    $length_insert = length ($insert);
                    calculate_frame_shift($length_insert, $insert);
            } elsif ($seq =~ /$id_lamda/) {
					$id = "lamda";
                	$end = $conserved_lamda;
                	$seq =~ /$start(.+)$end/;
                	$insert = $1;
	                $length_insert = length ($insert);
                	calculate_frame_shift($length_insert, $insert);
            } else {
                    $id = $fail;
                    $insert = $fail;
                    $length_insert = $fail;
                    $result = $fail;
                    $termination = $fail;
            }
    } elsif ($seq =~ /$global_revcomp/){
            revcompDNA($seq);
            print "reverse it\n";
    } else { 
			$leader = $fail;
			$id = $fail;
            $insert = $fail;
            $length_insert = $fail;
            $result = $fail;
            $termination = $fail;
    }
}

sub calculate_frame_shift{
    my ($value, $insert)= @_;
    if (defined($value)){
    $result= $value % 3 ;
        if ($result  == 0) {
            early_termination_check($insert);
        } else {
            $result = $fail;
        }
        push (@frame_shift, $result);
        } else { 
            $result = $fail;
            $length_insert = $fail;
        }
}

sub early_termination_check {

    my ($string) = @_;
#my $termination;
    # 'true' means early sequence terminates prematurely, 'false' means no stop codon found
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
    # looks for index of stop codon. Does modulus function to see if it is in the reading frame. If it is an early termination
    # sequence, returns 'true'
    my $result = index($string, $char, $offset);
    $check = $result % 3;

    if ($check == 0) {
            $termination = "TRUE";
            return;
    } else {
        while ($result != -1) {
            $check = $result % 3;
            if ($check == 0) {
           	  $termination = "TRUE";
              return;
            }
            $offset = $result + 1;
            $result = index($string, $char, $offset);
        }
        if (! defined($stop[$file_count])) {
            my $result_2 = index($string, $char_2, $offset_2);
            $check_2 = $result_2 % 3;
            if ($check_2 == 0) {
                    $termination = "TRUE";
                	return;
            } else {
                while ($result_2 != -1) {
                    $check_2 = $result_2 % 3;
                    if ($check_2 == 0) {
           				$termination = "TRUE";
                        return;
                    }
                    $offset_2 = $result_2 +1;
                    $result_2 = index($string, $char_2, $offset_2);
                }
            }
            if (! defined($stop[$file_count])){
                my $result_3 = index($string, $char_3, $offset_3);
                $check_3 = $result % 3;
                if ($check_3 == 0) {
            		$termination = "TRUE";
                    return;
                } else {
                    while ($result_3 != -1) {
                        $check_3 = $result_3 % 3;
                        if ($check_3 == 0) {
                    	    $termination = "TRUE";
                            return;
                        }
                        $offset_3 = $result_3 +1;
                        $result_3 = index($string, $char_3, $offset_3);
                    }
                }
            }
        }
        if (! defined($stop[$file_count])) {
          $termination = "FALSE";           
        }
    }
}

sub revcompDNA {
    my ($DNA) = @_;
    my $rev = reverse ($DNA); 
    $rev =~ tr/ACGTacgt/TGCAtgca/; 
    insert_id($rev, $const_frame_leader_all, $const_frame_leader_all_revcomp, $constant_vector_gamma, $constant_vector_kappa, $constant_vector_lamda, $const_frame_constant_gamma, $const_frame_constant_kappa, $const_frame_constant_lamda);	
}