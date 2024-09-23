#!/usr/bin/perl

#Sergio Vargas R. 2015. Munich, Germany

use strict;
use warnings;
use HTTP::Tiny;
use LWP::UserAgent;
use Getopt::Long;

#ADD GNU license here!

#command line options, the variable names should be self explanatory;
my $input_file;
my $verbose = "";
my $debug = "";
my $output_file = "out.tsv"; #where should I store the genomes?
my $save_tmps = "";
my $list_file = undef;
my $with_protein_annotations = "";


my $help = "Get_GO_Annotations.pl

The following options have to be provided:
		--in in.tsv; tab separated file to be read
		--out out.tsv; suffix for the output files?

these options can be passed to the program:
		--list list.txt; list entries to retrieve if the annotations of a full transcriptome are given as input
		--with_protein_annotations T; if you want to get the name of the annotated protein and of the taxon of the annotated protein (e.g. Mus musculus)
		--save_tmps F; if you want to keep the annotations on a gene level use T here.
 		--verbose F
 		--debug F
 		
A normal command line would look like:
	perl Get_GO_Annotations.pl -in ./file.tsv --list list.txt

The option --list requires you to pass a list of transcripts of interest which GOs will be fetched.

The list is simply a list of transcripts with no headers, footers or blank lines, e.g.:

TR4|c9_g1_i1
TR41535|c9_g1_i1
TR41536|c9_g1_i1

Note:

this script expects trinity transcripts annotated using blast agains SWISSPROT (with output 6). This looks like:

TR41535|c9_g1_i1	sp|P70414|NAC1_MOUSE	53.36	937	377	12	402	3065	45	970	0.0	  793

it extracts the SWISSPROT accesions (i.e. P70414) and looks for the GO terms associated with this accession.

The script produces three output files in the same folder where the script is:

function_out.tsv
process_out.tsv
component_out.tsv

with the go terms for each transcript as follows:

P70414	TR41535|c9_g1_i1	GO:XXXX,GO:YYYY

these files should be then usable for enrichment analyses pipelines like top_go.

If the flag --with_protein_annotations is active, the output files look like:

P70414	Protein_name_here	Taxon_name_here	TR41535|c9_g1_i1	GO:XXXX,GO:YYYY

";

#get options from command line
GetOptions("in=s" => \$input_file,
					 "out=s" => \$output_file,
					 "list=s" => \$list_file,
					 "with_protein_annotations" => \$with_protein_annotations,
					 "save_tmps=s" => \$save_tmps,
					 "verbose=s" => \$verbose,
					 "debug=s" => \$debug,
					 "help" => \&help
					 ) or die("Error in the command line arguments\n$help");

sub help{
	print $help;
}

sub Get_GO_Terms{

	my $ref_to_agent = shift;
	my $accession_to_search = shift;
	
	print "searching for ", $accession_to_search, "\n\n" if($verbose);

	#string to store GO terms for every "aspect", goIDs are separated by ","
	my $aspect_function = "";
	my $aspect_process = "";
	my $aspect_component = "";
	my $prot_name = "";
	my $tax_name = "";

	#string to returns one line with the "aspect" data for a gene of interest. Aspects are separated by #
	my $go_terms = "";

	#my $request_url = "http://www.ebi.ac.uk/QuickGO/GAnnotation?protein=$accession_to_search&format=tsv&col=proteinID,proteinSymbol,evidence,goID,goName,aspect,ref,with,from,proteinName,proteinTaxonName";

	my $request_url = "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?geneProductId=$accession_to_search";
	
	print $request_url, "\n" if $debug;

	print $accession_to_search, "\n" if $debug;
	
    my $tmp_filename = $accession_to_search . ".tmp";
	
	print $tmp_filename, "\n\n" if $debug;
	
	my $request = HTTP::Request->new(GET => $request_url, [ 'Accept' => 'text/tsv', 'Content_Disposition' => "form-data; name=attachment; filename=$tmp_filename"]);

	my $response = ${$ref_to_agent}->request($request);
    
	#note: this is not as clean as I would have like it to be but we are kind of in a rush.
	if($response->is_success){# I am assuming error will yield a 0 response...
        
		#code hacked from here https://www.ebi.ac.uk/QuickGO/WebServices.html
		open TMP, '>', $tmp_filename or die "Cannot open file $tmp_filename:", $!;
		
		print TMP $response->decoded_content;
		
		close TMP;
		
		print $response->as_string, "\n" if($debug);
		
		open TMP, '<', $tmp_filename or die "Cannot open file $tmp_filename:", $!;
		
		my $head = <TMP>; #discard the header.
		while (<TMP>) {
		
            #GENE_PRODUCT_DB|GENE_PRODUCT_ID|SYMBOL|QUALIFIER|GO_TERM|GO ASPECT|ECO_ID|GO_EVIDENCE|CODE|REFERENCE|WITH/FROM|TAXON_ID|ASSIGNED_BY|ANNOTATION_EXTENSION|DATE

			chomp;
			
			my ($gene_product, $protein_id, $proteinSymbol, $qualifier, $goID, $aspect, $eco_id, $evidence, $ref, $with_from, $protein_taxon_name) = split(/\t/);
			
			print $protein_id, "\t", $goID, "\t", $aspect, "\n" if $verbose;
			
			if($aspect eq "F"){

				$aspect_function .= "$goID,";

			}elsif($aspect eq "P"){

				$aspect_process .= "$goID,";

			}elsif($aspect eq "C"){

				$aspect_component .= "$goID,";
				
			}else{
				#for some reason, none of the above!
			}
			
			$prot_name = $protein_id;
			$tax_name = $protein_taxon_name;
		}
		close TMP;
		print "cleaning temporary files\n" if $verbose;
		unlink "$accession_to_search.tmp" or warn "cannot remove file $accession_to_search.tmp: $!" unless $save_tmps;
	}else{
        print "no answer...\n"; die;
	}
	
	print $aspect_function, "\n", $aspect_process, "\n", $aspect_component, "\n" if $verbose;
	$go_terms = $aspect_function . "#" . $aspect_process . "#" . $aspect_component;
	
	$go_terms = $prot_name . "#" . $tax_name . "#" . $go_terms if $with_protein_annotations;
	
	return $go_terms;
}

sub Hot_Print_Outfile{

	#selected passed filehandle for hot print and returns a print status.
	#before returning selects STDOUT back.

	my $filehandle = shift;
	my $string = shift;
	
	select $filehandle;
	$| = 1; #disable buffering on selected filehandle

	my $print_status = print $filehandle $string;
	
	$| = 0;#enable buffering
	
	select STDOUT;
	return $print_status;
}

sub Get_List_of_Queries{
	my $list_filehandle = shift;

	my @target_transcripts;
	my $targets=0;
	
	while(<$list_filehandle>){
		chomp;
		s/\|//g;
#		my ($id, $transcript) = split("\t");
#		print $transcript, "\n" if($debug);
#		push(@target_transcripts, $transcript);

		print $_, "\n" if($debug);
		push(@target_transcripts, $_);

		$targets++;
		}

	print "$targets in list file\n" if($verbose);

	return @target_transcripts;
}

sub main{

	#check that some required options are defined;
	die "\n$help\n" unless defined $input_file;
	
	print "opening $input_file\n" if $verbose;
	open INFILE, $input_file or die "Cannot open file $input_file:", $!;

	my @query_list;
	my $list_file_filehandle;
	
	if(defined $list_file){
		print "opening $list_file\n" if $verbose && $list_file;
		open $list_file_filehandle, '<' ,$list_file or die "Cannot open list file: ", $!;
		
		@query_list=Get_List_of_Queries($list_file_filehandle);
	}

	print "opening output files\n" if $verbose;
	open my $function_filehandle, '>', "function_$output_file" or die "Cannot open file for output: ", $!;
	open my $process_filehandle, '>', "process_$output_file" or die "Cannot open file for output: ", $!;
	open my $component_filehandle, '>', "component_$output_file" or die "Cannot open file for output: ", $!;

  open my $log_filehandle, '>', "Get_GO_Annotations.log" or die "Cannot open file for logging: ", $!;

	#create a new user agent and pass a ref to it to be used in sub routine
	my $ua = LWP::UserAgent->new;
	#my $ua = HTTP::Tiny->new();
	
	my $query_counter = 1;
	
	while(<INFILE>){
		print if $debug;
		
		chomp;
		#s/\t/\|/g;

		print if $debug;

		my @fields = split(/\t/);#SWISSPROT accession number is in field 1 in the array

		print "\n", $fields[0], $fields[1], "\n" if($debug);
		
		#my $transcript_id = join("",$fields[0], $fields[1]);
		my $transcript_id = $fields[0];
		my $uniprot_accession = $fields[1];
		$transcript_id =~ s/_i[0-9]*.p[0-9]*//g;
        $uniprot_accession =~ s/\.[0-9]*$//g;
		
		print "\n", $transcript_id, "\t", $uniprot_accession, "\n" if($verbose);
		
    my $in_list = 1;

    if(defined $list_file){
      $in_list = grep $_ eq $transcript_id, @query_list;
    }
		
		if($in_list){
			print "$transcript_id in query list: $in_list\nAbout to fetch GO terms\n" if($verbose);
		
			print STDERR "\n$transcript_id: Querying SWISSPROT accession $uniprot_accession...Total fetched $query_counter\r";
		
		  my $log_line = "$transcript_id fetched\n";
		  warn "Error printing to log file!" if(!&Hot_Print_Outfile($log_filehandle, $log_line));
		
			print "\n", $fields[3], "\n" if $debug;
		
			my ($protein, $taxon, $function_go_terms, $process_go_terms, $component_go_terms, $function_out_line, $process_out_line, $component_out_line) = "NA";
			
			if($with_protein_annotations){
			  ($protein, $taxon, $function_go_terms, $process_go_terms, $component_go_terms) = split(/#/, Get_GO_Terms(\$ua, $fields[3]));
			  
			  $function_out_line = $fields[3] . "\t" . $protein . "\t" . $taxon . "\t" . $fields[0] . "|" . $fields[1] . "\t" . substr($function_go_terms, 0, -1) . "\n";
		  	warn "Error printing GO Function output!" if(!&Hot_Print_Outfile($function_filehandle, $function_out_line));
			  #print FUNCTION $function_out_line or warn "writing function failed\n";

			  $process_out_line =  $fields[3] . "\t" . $protein . "\t" . $taxon . "\t" . $fields[0] . "|" . $fields[1] . "\t" . substr($process_go_terms, 0, -1) . "\n";
			  warn "Error printing GO Process output!" if(!&Hot_Print_Outfile($process_filehandle, $process_out_line));
		  	#print PROCESS $process_out_line or warn "writing process failed\n";

			  $component_out_line = $fields[3] . "\t" . $protein . "\t" . $taxon . "\t" . $fields[0] . "|" . $fields[1] . "\t" . substr($component_go_terms, 0, -1) . "\n";
			  warn "Error printing GO Component output!" if(!&Hot_Print_Outfile($component_filehandle, $component_out_line));
  			#print COMPONENT $component_out_line or warn "writing component failed\n";			  
			}else{
			  ($function_go_terms, $process_go_terms, $component_go_terms) = split(/#/, Get_GO_Terms(\$ua, $uniprot_accession));
			  
			  $function_out_line = $uniprot_accession . "\t" . $fields[0] . "\t" . substr($function_go_terms, 0, -1) . "\n";
		  	warn "Error printing GO Function output!" if(!&Hot_Print_Outfile($function_filehandle, $function_out_line));
			  #print FUNCTION $function_out_line or warn "writing function failed\n";

			  $process_out_line =  $uniprot_accession . "\t" . $fields[0] . "\t" . substr($process_go_terms, 0, -1) . "\n";
			  warn "Error printing GO Process output!" if(!&Hot_Print_Outfile($process_filehandle, $process_out_line));
		  	#print PROCESS $process_out_line or warn "writing process failed\n";

			  $component_out_line = $uniprot_accession . "\t" . $fields[0] . "\t" . substr($component_go_terms, 0, -1) . "\n";
			  warn "Error printing GO Component output!" if(!&Hot_Print_Outfile($component_filehandle, $component_out_line));
  			#print COMPONENT $component_out_line or warn "writing component failed\n";			  
			}
			

			$query_counter++;
		}else{
			print "$transcript_id NOT in query list: $in_list\n" if $verbose;
		}
	}
	close INFILE;
	close $function_filehandle;
	close $process_filehandle;
	close $component_filehandle;
	close $log_filehandle;
	return 0;
}

exit &main();
