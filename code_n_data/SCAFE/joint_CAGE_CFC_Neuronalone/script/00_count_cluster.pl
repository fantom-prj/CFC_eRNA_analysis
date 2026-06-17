#!/usr/bin/perl -w
use strict;

my $baseDir = $ARGV[0];
my $ctss_scope_bed_path = $ARGV[1];
my $in_ctss_list = $ARGV[2];
my $countRegion_bed_path = $ARGV[3];
my $tag = $ARGV[4];

my $outDir = "$baseDir/per_lib_cluster/";
my $matrixDir = "$baseDir/count_matrix_cluster/";
my $run_scafe_script_path = '/osc-fs_home/hon-chun/analysis/tenX_single_cell/scafe/dev/deploy/release/v1.0.1/scripts/scafe.tool.bk.count';
my $genome = "hg38.gencode_v39";

system "mkdir -pm 755 $outDir";
system "mkdir -pm 755 $matrixDir";

my $lib_info_hsh_ref = {};
open (CTSSLIST, "<", $in_ctss_list);
while (<CTSSLIST>) {
	chomp;
	my ($libID, $ctss_bed_path) = split /\t/;
	die "ctss_bed_path $ctss_bed_path does not exists " if not -s $ctss_bed_path;
	$lib_info_hsh_ref->{$libID}{'ctss_bed_path'} = $ctss_bed_path;
}
close CTSSLIST;

my $data_hsh_ref = {};
my @lib_ary = sort keys %{$lib_info_hsh_ref};
foreach my $libID (@lib_ary) {
	my $ctss_bed_path = $lib_info_hsh_ref->{$libID}{'ctss_bed_path'};
	my $outputPrefix = $libID;
	my $force_rerun = 'no';
	my $max_thread = 10;
	my $count_log_path = "$outDir/$outputPrefix/log/count.log.tsv";
	$lib_info_hsh_ref->{$libID}{'count_log_path'} = $count_log_path;
	
	my $cmd = join " ", (
		"$run_scafe_script_path",
		"--countRegion_bed_path=$countRegion_bed_path",
		"--ctss_bed_path=$ctss_bed_path",
		"--ctss_scope_bed_path=$ctss_scope_bed_path",
		"--ctss_scope_slop_bp=0",
		"--genome=$genome",
		"--outputPrefix=$outputPrefix",
		"--outDir=$outDir",
	);
	system "$cmd" if not -s $count_log_path;
	
	open COUNTLOG, "<", $count_log_path;
	<COUNTLOG>;
	while (<COUNTLOG>) {
		chomp;
		my ($CREID, $count) = split /\t/;
		$data_hsh_ref->{$CREID}{$libID} = $count;
	}
	close COUNTLOG;
}

open OUTTABLE, ">", "$matrixDir/$tag.count.txt";
print OUTTABLE join "", (join "\t", ('CREID', @lib_ary)), "\n";
foreach my $CREID (sort keys %{$data_hsh_ref}) {
	my @output_ary = ($CREID);
	foreach my $libID (@lib_ary) {
		push @output_ary, $data_hsh_ref->{$CREID}{$libID};
	}
	print OUTTABLE join "", (join "\t", (@output_ary)), "\n";
}
close OUTTABLE;
