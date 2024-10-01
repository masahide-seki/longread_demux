#!/usr/bin/perl
#$ -S /usr/bin/perl
# Revision: 1.0.0
#  Date: 2024/03/29
#  Copyright c 2024, by DYNACOM Co.,Ltd.

#############################################################################################


use strict;
use Getopt::Long qw(:config posix_default no_ignore_case gnu_compat);
use FileHandle;
use Data::Dumper;

################################################################################
### Usage
#################################################################################
my %options=("distance" => 0 );
my @optary = ( "help|h" , "out|o:s", "input|i:s","list|l:s", "distance|d:i");
GetOptions(\%options, @optary) or &usage(0);

if ( $options{help} ) { &usage(1); }
sub usage {
	my ($flag, $mes) = @_;
	chomp(my $program = `basename $0`);
	my $usage = " Usage: $program -i input.fq.gz -l index.csv -o output  -d 100 \n";
	$usage   .= "        qsub -cwd $program -i input.fq.gz -l index.csv -o output -d 100 \n";
	$usage   .= "   e.g. $program \n \n";

	my $help = <<'&EOT&';
	- - - - - - - -
		[Options]
		-h --help  ... help
		-o --out   ... output file prefix            e.g sample
		-i --input   ... input fastq file .gz or not
		-l --list  ... barcode index list     
		-d --distance ... index position limit from read end

判定の流れ
・i7 / i5 いずれかが複数種類のID検出→Ambiguous
・i7 / i5 両方1種類のID 検出 ( 1種類が複数回検出も含む)
        →i7/i5の組み合わせが非一致→Other
        →i7/i5の向きが違う(fw,fw/rv,rv/ )→sample
　　　　→i7/i5が組み合わせ向きともに一致(fw/rv & rv/fw)　→sample
・i7 / i5 片方のみ1種類のID検出 ( 1種類が複数回検出も含む)
        →i7/i5が1種類　→sample 該当→ sample
・検出無し→NotDetect

&EOT&
	print " !! $mes\n\n" if $mes;
	print $usage;
	print $help if $flag;
	exit 1;
}

################################################################################
## Sub
#################################################################################

sub GetIndextStats{
	my $hashref = shift;
	my %hash = %{$hashref};

	my @detcted_list;
	my $num_of_detectIndex = 0;
	foreach my $index_name (keys( %hash )){
		$num_of_detectIndex++;
		foreach my $strand ("Fw","Rv" ){
			foreach my $pos (keys( %{$hash{$index_name}{$strand}} )){
				my $dist = $hash{$index_name}{$strand}{$pos};
				push( @detcted_list ,"$index_name,$strand,$pos,$dist" );
			}
		}
	}
	return( $num_of_detectIndex , join(":",@detcted_list) );
}

################################################################################
## main funciton
#################################################################################

sub main {
	my $options = shift;
	my $fq    = $$options{input};
	my $out   = $$options{out};
	my $index = $$options{list};
	my $dist  = $$options{distance};
	#print "$dist\n";
	
	#############################################
	#--- Get Sample Index infomation % make output file handle
	#############################################
	my %index_dt ;
	my %out_fq_handle ;
	my @summary_head ;
	open( my $s_infh , "< $index" ) || die "cannot open $index :$!";
		
	while (my $line =<$s_infh>){
		next if( $line =~/^#/);
		$line =~ s/\r\n|\n|\r//g;
		my ($sample ,$i7, $i5,$i7fw,$i7rv,$i5fw,$i5rv ) = split(/\t/,$line);
		$index_dt{i7}{$i7}{Fw} = $i7fw ;
		$index_dt{i5}{$i5}{Fw} = $i5fw ;
		$index_dt{i7}{$i7}{Rv} = $i7rv ;
		$index_dt{i5}{$i5}{Rv} = $i5rv ;
		$index_dt{i7}{$i7}{Sample} = $sample ;
		$index_dt{i5}{$i5}{Sample} = $sample ;
		open (  $out_fq_handle{$sample}   ,"> Demulti_${i7}-${i5}_${sample}.fq") || die "canot open Demulti_${i7}-${i5}_${sample}.fq :$!";
		push(@summary_head, $sample);
	}

	# other index detection
	open ($out_fq_handle{"NotDetect"}  ,"> Demulti_NotDetect.fq")  || die "canot open Demulti_NotDetect.fq  :$!";
	open ($out_fq_handle{"OtherIndex"} ,"> Demulti_OtherIndex.fq") || die "canot open Demulti_OtherIndex.fq :$!";
	open ($out_fq_handle{"Ambiguous"} ,"> Demulti_Ambiguous.fq")  || die "canot open Demulti_Ambiguous.fq  :$!";
	push(@summary_head, ( "OtherIndex","Ambiguous" ,"NotDetect"));
	close($s_infh );
	open ( my $oh1  ,"> logBcDetected_${out}_detail.txt") || die "canot open ${out}_R1.fq :$!";
	open ( my $oh2  ,"> logBcDetected_${out}_summary.txt") || die "canot open ${out}_R1.fq :$!";
	print $oh1 "##options: -i $fq -l $index -d $dist -o $out  ....  \n";
	print $oh1 join("\t","id","Detected_Index","Sample","length","Detect_i7","Detect_i5","Skip_i7","Skip_i7" ),"\n";

	#############################################
	#--- read input file
	#############################################
	my $fh;
	if ( $fq =~/.gz$/  ) {
		open ( $fh ,"zcat $fq| ") || die "canot open $fq :$!";
	} else  {
		open ( $fh ,"zcat $fq| ") || die "canot open $fq :$!";
	}
	my %summary;

	while (my $line = <$fh>){
		#last if( $count > 5000);
		$line=~s/\n|\r|\n\r//g;
		(my $line2 = <$fh>)=~s/\n|\r|\n\r//g ;
		(my $line3 = <$fh>)=~s/\n|\r|\n\r//g ;
		(my $line4 = <$fh>)=~s/\n|\r|\n\r//g ;
		my $read_id = (split(/\s/,$line))[0];
		$read_id =~s/^@//g;
		my %idx_jduge ;
		my $readLength =length($line2);

		# Index seq searche
		my @idxP_list = ("i7","i5"); 
		foreach my $idxP( @idxP_list ) {
			my @index_ids = keys( %{$index_dt{$idxP}});
			#print "@index_ids\n";
			foreach my $name( @index_ids ) {
				foreach my $strand ( ("Fw","Rv") ) {
					my $bc_seq = $index_dt{$idxP}{$name}{$strand};
					if( $line2 =~/$bc_seq/){
						my $bc_locus = -1;
						while(1){
							$bc_locus = index( $line2 , $bc_seq , $bc_locus+ 1);
							last if $bc_locus == -1;
							my $dist_from_Readend = ( $bc_locus <  ( $readLength - $bc_locus) ) ? $bc_locus :  ( $readLength - $bc_locus);
							my $dist_dlg =  ($dist != 0 & $dist_from_Readend > $dist  ) ? "SKIP":"OK";
							$idx_jduge{$idxP}{$dist_dlg}{$name}{$strand}{$bc_locus} = $dist_from_Readend;
						}
					}
				}
			}
		}
		#############################################
		# judge
		#############################################
		my $flg = "NotDetect";
		my $detect_index = "NotDetect";
		my $Sample_Naems = "-";

		my ($num_of_detectI7,$Detect_i7) = (exists( $idx_jduge{"i7"}{OK})) ? GetIndextStats( $idx_jduge{"i7"}{OK} ) : (0,"-") ;
		my ($num_of_detectI5,$Detect_i5) = (exists( $idx_jduge{"i5"}{OK})) ? GetIndextStats( $idx_jduge{"i5"}{OK} ) : (0,"-") ;
		my (undef,$Skip_i7) = (exists( $idx_jduge{"i7"}{SKIP})) ? GetIndextStats( $idx_jduge{"i7"}{SKIP} ) : (0,"-") ;
		my (undef,$Skip_i5) = (exists( $idx_jduge{"i5"}{SKIP})) ? GetIndextStats( $idx_jduge{"i5"}{SKIP} ) : (0,"-") ;

		my @i7_detected = ( $Detect_i7 eq "-" ) ? () :split(/:/ , $Detect_i7 );
		my @i5_detected = ( $Detect_i5 eq "-" ) ? () :split(/:/ , $Detect_i5 );

		if ( $num_of_detectI7 > 1 or $num_of_detectI5 > 1) {
			$flg = "Ambiguous";
			$detect_index ="Ambiguous";
		} elsif ( $num_of_detectI7 == 1 and $num_of_detectI5 == 1 ){ #i7 & i5
			# sample name
			my ( $i7name , $i7st , $i7pos , $i7dist ) = split(/,/ ,$i7_detected[0]);
			my ( $i5name , $i5st , $i5pos , $i5dist ) = split(/,/ ,$i5_detected[0]);
			my $i7_sampleID = $index_dt{"i7"}{$i7name}{Sample} ;
			my $i5_sampleID = $index_dt{"i5"}{$i5name}{Sample} ;

			$flg = "OtherIndex";
			$detect_index = "$i7name-$i5name";
			if ( ( $i7_sampleID ne "") and (  $i7_sampleID eq $i5_sampleID ) ){
				$flg = $i7_sampleID;
				$Sample_Naems = $i7_sampleID;
			}
	
		} elsif( $num_of_detectI7 == 1 or $num_of_detectI5 == 1 ) { #i7 
			# my @detect_index_list = keys( $idx_jduge{"i7"} );
			my $use_inxP = ( $num_of_detectI7 == 1 ) ? $i7_detected[0]: $i5_detected[0];
			my ( $iXname , $iXst , $iXpos , $iXdist ) = split(/,/ ,$use_inxP);
			my $iX_sampleID = ( $num_of_detectI7 == 1 ) ? $index_dt{"i7"}{$iXname}{Sample}:$index_dt{"i5"}{$iXname}{Sample};

			$flg = "OtherIndex";
			$detect_index = "$iXname";
			if( $iX_sampleID ne "" ){
				$flg = $iX_sampleID;
				$Sample_Naems = $iX_sampleID;
			}
		}
		$summary{$flg}++;

		my $out_fh  = $out_fq_handle{$flg};
		print $oh1 join("\t",$read_id , $detect_index,$Sample_Naems, $readLength ,$Detect_i7,$Detect_i5,$Skip_i7,$Skip_i5 ),"\n";
		print $out_fh "\@$read_id $detect_index:$Sample_Naems\n$line2\n$line3\n$line4\n";
	}

	print $oh2 join("\t", @summary_head),"\n";
	my @sum_val =map{ $summary{$_}} @summary_head;
	print $oh2 join("\t", @sum_val),"\n";

	close($fh);
	close($oh1);close($oh2);#close($ohfq);
	foreach my $sample (keys( %out_fq_handle )){
		close( $out_fq_handle{$sample}  );
	}
}

################################################################################
## Main
#################################################################################
main( \%options);


exit;

