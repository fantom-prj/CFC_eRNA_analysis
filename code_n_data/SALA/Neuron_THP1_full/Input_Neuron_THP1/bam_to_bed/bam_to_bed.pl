#!/usr/bin/perl -w
use strict;

my $outDir = "/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/SCAFE/20231110/Neuron_THP1/ontCAGE/bam_to_bed";
my $in_lib_list_path = "/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/data_link/dorado.Neuron_THP1.BAM.txt";
my $genome = "hg38.gencode_v39";
my $dirPath = "/osc-fs_home/hon-chun/analysis/tenX_single_cell/scafe/dev/deploy/release/v1.0.1/scripts/";

	my $tabix_bin = "$dirPath/../resources/bin/tabix/tabix";
	my $bgzip_bin = "$dirPath/../resources/bin/bgzip/bgzip";

my $bed_dir = "$outDir/bed/";
my $subsample_100x_bed_dir = "$outDir/subsample_100x_bed/";
my $subsample_10x_bed_dir = "$outDir/subsample_10x_bed/";
my $sbatch_dir = "$outDir/00_sbatch/";
system "mkdir -pm 755 $outDir";
system "mkdir -pm 755 $sbatch_dir";
system "mkdir -pm 755 $bed_dir";
system "mkdir -pm 755 $subsample_100x_bed_dir";
system "mkdir -pm 755 $subsample_10x_bed_dir";
system "cp $0 $outDir";

my $start_time_log_path = "$outDir/00_start.time.log.txt";
my $finish_time_log_path = "$outDir/00_finish.time.log.txt";
system"echo \"========== All runs are started at \$(date)\" ========== >$start_time_log_path\n";
system"echo \"========== All runs are started at \$(date)\" ========== >$finish_time_log_path\n";

my $lib_info_hsh_ref = {};

open BAMLIST, "<", $in_lib_list_path;
while (<BAMLIST>) {
	chomp;
	next if $_ =~ m/^#/;
	my ($libID, $prefix, $bamPath) = split /\t/;
	
	die "bamPath does not exists for $libID\n" if not -s $bamPath;
	my $bedPath = "$bed_dir/$libID.bed.bgz";
	my $subsample_100x_bedPath = "$subsample_100x_bed_dir/$libID.bed.bgz";
	my $subsample_10x_bedPath = "$subsample_10x_bed_dir/$libID.bed.bgz";
	$lib_info_hsh_ref->{$libID}{'prefix'} = $prefix;
	$lib_info_hsh_ref->{$libID}{'bedPath'} = $bedPath;
	$lib_info_hsh_ref->{$libID}{'subsample_100x_bedPath'} = $subsample_100x_bedPath;
	$lib_info_hsh_ref->{$libID}{'subsample_10x_bedPath'} = $subsample_10x_bedPath;
	$lib_info_hsh_ref->{$libID}{'bamPath'} = $bamPath;
}
close BAMLIST;

my $num_lib = keys %{$lib_info_hsh_ref};
print "$num_lib libraries read\n";

foreach my $libID (sort keys %{$lib_info_hsh_ref}) {
	my $bamPath = $lib_info_hsh_ref->{$libID}{'bamPath'};
	my $bedPath = $lib_info_hsh_ref->{$libID}{'bedPath'};
	my $prefix = $lib_info_hsh_ref->{$libID}{'prefix'};
	my $subsample_100x_bedPath = $lib_info_hsh_ref->{$libID}{'subsample_100x_bedPath'};
	my $subsample_10x_bedPath = $lib_info_hsh_ref->{$libID}{'subsample_10x_bedPath'};
		
	my $sbatch_sh_path = "$sbatch_dir/$libID.sbatch.cmd.sh";
	my $sbatch_stderr_path = "$sbatch_dir/$libID.sbatch.stderr.txt";
	my $sbatch_stdout_path = "$sbatch_dir/$libID.sbatch.stdout.txt";

	open (SBATCHCMD, ">", $sbatch_sh_path);
	print SBATCHCMD "#!/bin/bash\n";
	print SBATCHCMD "echo \"$libID is started at \$(date)\" >>$start_time_log_path\n";
	print SBATCHCMD "samtools view -bh $bamPath | bedtools bamtobed -i stdin -bed12 | sed 's/$libID/$prefix/g' | sort --parallel 10 -k1,1 -k2,2n | $bgzip_bin --threads 10 -c >$bedPath\n";
	print SBATCHCMD "$tabix_bin -p bed $bedPath\n";
	print SBATCHCMD "$bgzip_bin -dc $bedPath | awk 'NR % 100 == 0' | $bgzip_bin --threads 10 -c >$subsample_100x_bedPath\n";
	print SBATCHCMD "$tabix_bin -p bed $subsample_100x_bedPath\n";
	print SBATCHCMD "$bgzip_bin -dc $bedPath | awk 'NR % 10 == 0' | $bgzip_bin --threads 10 -c >$subsample_10x_bedPath\n";
	print SBATCHCMD "$tabix_bin -p bed $subsample_10x_bedPath\n";
	print SBATCHCMD "echo \"$libID is finished at \$(date)\" >>$finish_time_log_path\n";
	close SBATCHCMD;

	print "$libID sbatch submitted.\n";
	system "sbatch -c 5 -e $sbatch_stderr_path -o $sbatch_stdout_path $sbatch_sh_path";
}
