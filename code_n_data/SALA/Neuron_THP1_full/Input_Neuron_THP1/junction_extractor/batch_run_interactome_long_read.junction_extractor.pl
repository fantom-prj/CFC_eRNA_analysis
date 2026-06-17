#!/usr/bin/perl -w
use strict;

my $outDir = "/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/Input_Neuron_THP1/junction_extractor";
my $in_lib_list_path = "/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/data_link/dorado.Neuron_THP1.b4.transcriptClean.BAM.txt";
my $tabix_bin='/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SCAFE/resources/bin/tabix/tabix';
my $bgzip_bin='/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SCAFE/resources/bin/bgzip/bgzip';
my $samtools_bin='/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SCAFE/resources/bin/samtools/samtools';
my $bedtools_bin='/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SCAFE/resources/bin/bedtools/bedtools';
my $chrom_size_path='/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SCAFE/resources/genome/hg38.gencode_v39/tsv/chrom.sizes.tsv';
my $chrom_fasta_path='/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SCAFE/resources/genome/hg38.gencode_v39/fasta/genome.fa';
my $junction_extractor_script = '/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/perl_script_for_SALA/junction_extractor/junction_extractor_v0.1.pl';
my $max_thread = 20;
my $min_nt_qual = 10;
my $min_MAPQ = 20;
my $sbatch_dir = "$outDir/00_sbatch/";
system "mkdir -pm 755 $outDir";
system "mkdir -pm 755 $sbatch_dir";
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
	my ($libID, $bamPath) = split /\t/;
	
	die "bamPath does not exists for $libID\n" if not -s $bamPath;
	$lib_info_hsh_ref->{$libID}{'bamPath'} = $bamPath;
}
close BAMLIST;

my $num_lib = keys %{$lib_info_hsh_ref};
print "$num_lib libraries read\n";

foreach my $libID (sort keys %{$lib_info_hsh_ref}) {
	my $in_bam = $lib_info_hsh_ref->{$libID}{'bamPath'};
	my $out_prefix = $libID;
	my $out_dir = "$outDir/output/";
	my $sbatch_sh_path = "$sbatch_dir/$libID.sbatch.cmd.sh";
	my $sbatch_stderr_path = "$sbatch_dir/$libID.sbatch.stderr.txt";
	my $sbatch_stdout_path = "$sbatch_dir/$libID.sbatch.stdout.txt";

	my $junction_extractor_cmd = join " ", (
		"perl $junction_extractor_script",
		"--in_bam=".$in_bam,
		"--chrom_size_path=".$chrom_size_path,
		"--chrom_fasta_path=".$chrom_fasta_path,
		"--out_prefix=".$out_prefix,
		"--out_dir=".$out_dir,
		"--max_thread=".$max_thread,
		"--min_nt_qual=".$min_nt_qual,
		"--min_MAPQ=".$min_MAPQ,
		"--samtools_bin=".$samtools_bin,
		"--bedtools_bin=".$bedtools_bin,
		"--tabix_bin=".$tabix_bin,
		"--bgzip_bin=".$bgzip_bin,
	);

	open (SBATCHCMD, ">", $sbatch_sh_path);
	print SBATCHCMD "#!/bin/bash\n";
	print SBATCHCMD "echo \"$libID is started at \$(date)\" >>$start_time_log_path\n";
	print SBATCHCMD "$junction_extractor_cmd\n";
	print SBATCHCMD "echo \"$libID is finished at \$(date)\" >>$finish_time_log_path\n";
	close SBATCHCMD;

	print "$libID sbatch submitted.\n";
	system "sbatch -c 10 -e $sbatch_stderr_path -o $sbatch_stdout_path $sbatch_sh_path";
}
