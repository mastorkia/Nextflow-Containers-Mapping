use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use Log::Log4perl qw(:easy);

use Cwd  qw(abs_path); 

#Usage: perl generate_adapter_statistics.pl RTS.Pool2.list
# List contains files in the format:/share/data/IlluminaData/RTS/Mar_2016/Pool1/Primer_trimmed/error_log/RTS-0038_AAGGAT_L001.Nextera_PCR_primer_rc.trimmed_log.o85980 


my $file_list = shift (@ARGV);
unless (defined $file_list) {
        print "1.Please provide the sample file list containing complete path names of cutadapt .o files\n";
	exit;
}

#my $logger = Log::Log4perl->get_logger('Starting cutadapt statistics generation for $file_list');
open (LOGS, "<$file_list");

print "Sample_Name\tAdaptor_Name\tTotal_Reads_Processed\tNum_Reads_with_adaptor\tPercent_Reads_with_adaptor\tCutadapt_Parameters\n";

my %sample_hash;
my $count = 0;

while(<LOGS>) {
	chomp;
	my $log_file = $_;
	my @names = split(/\./, $_);
	my $cutadapt_info = `head -n 2 $log_file`;
	my @cutadapt_info_array = split(/\n/, $cutadapt_info);
 	
	my $cutadapt_cmd = $cutadapt_info_array[1];
	my $cutadapt_params = "";
	
	if ($cutadapt_cmd =~ /.*\:\s(.*)\-\-.*/ ) {
		$cutadapt_params = $1;
	}
	my $sample_loc = $names[0];
	my ($sample, $dir) = fileparse($sample_loc);
	#print "$sample\t";
	my $adaptor = $names[1];
	
	my $read_counts = `grep \"Total reads processed:\" $log_file`;
	my @read_count_array = split(/\s{2,}/, $read_counts);
	chomp($read_count_array[1]);
	
	my $trimmed_counts = `grep \"Reads with adapters:\" $log_file`;
	my @trimmed_count_array = split(/\s{2,}/, $trimmed_counts);
	my $trimmed_stats = $trimmed_count_array[1];
	my @trimmed_no_percent = split(/\s+/, $trimmed_stats);
	
	my $num_reads_trimmed = $trimmed_no_percent[0];
	my $percent_reads_trimmed = $trimmed_no_percent[1];	
	$percent_reads_trimmed =~ s/\(//;
	$percent_reads_trimmed =~ s/\)//;	

	my $string_to_print = "$adaptor\t$read_count_array[1]\t$num_reads_trimmed\t$percent_reads_trimmed\t$cutadapt_params\n";
	#print "$string_to_print";
	push @{$sample_hash{$sample}}, $string_to_print;
	}

	foreach my $key (sort keys %sample_hash) {
		for my $i ( 0 .. $#{ $sample_hash{$key} } ) {
        		print "$key\t$sample_hash{$key}[$i]";
    		}
    		print "\n";
	} 



