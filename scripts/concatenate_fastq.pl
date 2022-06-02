use strict;
use strict;
use warnings;
use Data::Dumper;
use File::Basename;

my $file_list = shift @ARGV; # || die "Please provide fastq files for concatenation";

my %sample_name_hash;

open(LIST, "<$file_list");
my @array = ();
my $dir = "";

while(<LIST>) { 
	chomp;
	my $this_fastq = $_;
	(my $file, $dir) = fileparse($this_fastq);
	#print "$file\n$dir\n";
	my @names = split("_", $file);
	my $sample = join("_", @names[0..1]);
	#print $sample, "\n";
	push @{$sample_name_hash{$sample}}, $this_fastq;

}

print $dir, "\n";

#print Dumper %sample_name_hash;


foreach my $sample_name (keys %sample_name_hash) {
	print $sample_name, "\n";
	foreach my $fastq_file (@{$sample_name_hash{$sample_name}}) { 
		my $cmd = "cat $fastq_file >> $dir$sample_name.fastq";
		print $cmd, "\n";
		system("$cmd");
	}
}
