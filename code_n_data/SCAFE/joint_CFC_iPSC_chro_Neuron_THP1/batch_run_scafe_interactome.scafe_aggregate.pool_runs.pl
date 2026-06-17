#!/usr/bin/perl -w
use strict;

my $baseDir='/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/SCAFE/20231126/iPSchro/ontCAGE/aggregate/pool_runs';
my $run_script = '/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/SCAFE/20231126/iPSchro/ontCAGE/script/00_run_aggregate_bam_to_ctss_only.sh';
my $run_list_path = '/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/SCAFE/20231126/iPSchro/ontCAGE/aggregate/00_iPSchro.aggregate_list.pool_run.txt';
my $genome = "hg38.gencode_v39";
my $sbatch_dir = "$baseDir/00_sbatch/";
my $aggregate_list_dir = "$baseDir/00_aggregate_list/";
my $ctss_dir = "$baseDir/ctss_bed/";
my $all_ctss_dir = "$ctss_dir/all_ctss/";
my $ung_ctss_dir = "$ctss_dir/ung_ctss/";
my $run_info_hsh_ref = {};
system "mkdir -pm 755 $baseDir";
system "mkdir -pm 744 $sbatch_dir";
system "mkdir -pm 744 $aggregate_list_dir";
system "mkdir -pm 744 $all_ctss_dir";
system "mkdir -pm 744 $ung_ctss_dir";
system "cp $0 $baseDir";
open ALLPOOLCTSSLIST, ">", "$ctss_dir/00_pool_runs.ctss_list.tsv";

my $start_time_log_path = "$baseDir/00_start.time.log.txt";
my $finish_time_log_path = "$baseDir/00_finish.time.log.txt";
system"echo \"========== All runs are started at \$(date)\" ========== >$start_time_log_path\n";
system"echo \"========== All runs are started at \$(date)\" ========== >$finish_time_log_path\n";

open LIST, "<", $run_list_path;
while (<LIST>) {
	#ontCAGE.iPSC.rep1	ontCAGE.iPSC.rep1_run1	/osc-fs_home/hon-chun/analysis/FANTOM6/interactome/all_CAGE/scafe/20221226/ctss_bed/all_ctss//ontCAGE.iPSC.rep1_run1.all.ctss.bed.gz	/osc-fs_home/hon-chun/analysis/FANTOM6/interactome/all_CAGE/scafe/20221226/ctss_bed/ung_ctss//ontCAGE.iPSC.rep1_run1.ctss.bed.gz
	#ontCAGE.iPSC.rep1	ontCAGE.iPSC.rep1_run2	/osc-fs_home/hon-chun/analysis/FANTOM6/interactome/all_CAGE/scafe/20221226/ctss_bed/all_ctss//ontCAGE.iPSC.rep1_run2.all.ctss.bed.gz	/osc-fs_home/hon-chun/analysis/FANTOM6/interactome/all_CAGE/scafe/20221226/ctss_bed/ung_ctss//ontCAGE.iPSC.rep1_run2.ctss.bed.gz
	#ontCAGE.iPSC.rep1	ontCAGE.iPSC.rep1_run3	/osc-fs_home/hon-chun/analysis/FANTOM6/interactome/all_CAGE/scafe/20221226/ctss_bed/all_ctss//ontCAGE.iPSC.rep1_run3.all.ctss.bed.gz	/osc-fs_home/hon-chun/analysis/FANTOM6/interactome/all_CAGE/scafe/20221226/ctss_bed/ung_ctss//ontCAGE.iPSC.rep1_run3.ctss.bed.gz
	#ontCAGE.iPSC.rep2	ontCAGE.iPSC.rep2_run1	/osc-fs_home/hon-chun/analysis/FANTOM6/interactome/all_CAGE/scafe/20221226/ctss_bed/all_ctss//ontCAGE.iPSC.rep2_run1.all.ctss.bed.gz	/osc-fs_home/hon-chun/analysis/FANTOM6/interactome/all_CAGE/scafe/20221226/ctss_bed/ung_ctss//ontCAGE.iPSC.rep2_run1.ctss.bed.gz
	#ontCAGE.iPSC.rep2	ontCAGE.iPSC.rep2_run2	/osc-fs_home/hon-chun/analysis/FANTOM6/interactome/all_CAGE/scafe/20221226/ctss_bed/all_ctss//ontCAGE.iPSC.rep2_run2.all.ctss.bed.gz	/osc-fs_home/hon-chun/analysis/FANTOM6/interactome/all_CAGE/scafe/20221226/ctss_bed/ung_ctss//ontCAGE.iPSC.rep2_run2.ctss.bed.gz
	#ontCAGE.iPSC.rep2	ontCAGE.iPSC.rep2_run3	/osc-fs_home/hon-chun/analysis/FANTOM6/interactome/all_CAGE/scafe/20221226/ctss_bed/all_ctss//ontCAGE.iPSC.rep2_run3.all.ctss.bed.gz	/osc-fs_home/hon-chun/analysis/FANTOM6/interactome/all_CAGE/scafe/20221226/ctss_bed/ung_ctss//ontCAGE.iPSC.rep2_run3.ctss.bed.gz
	chomp;
	my ($runID, $libID, $all_ctss_path, $ung_ctss_path) = split /\t/;
	$run_info_hsh_ref->{$runID}{$libID}{'all_ctss_path'} = $all_ctss_path;
	$run_info_hsh_ref->{$runID}{$libID}{'ung_ctss_path'} = $ung_ctss_path;
}
close LIST;

foreach my $runID (sort keys %{$run_info_hsh_ref}) {
	my $lib_list_path = "$aggregate_list_dir/00_$runID.aggregate_list.tsv";
	open AGRLIST, ">", $lib_list_path;
	foreach my $libID (sort keys %{$run_info_hsh_ref->{$runID}}) {
		my $all_ctss_path = $run_info_hsh_ref->{$runID}{$libID}{'all_ctss_path'};
		my $ung_ctss_path = $run_info_hsh_ref->{$runID}{$libID}{'ung_ctss_path'};
		print AGRLIST join "", (join "\t", ($libID, $all_ctss_path, $ung_ctss_path)), "\n";
	}
	close AGRLIST;
	my $outputPrefix = $runID;
	my $sbatch_sh_path = "$sbatch_dir/$outputPrefix.sbatch.cmd.sh";
	my $sbatch_stderr_path = "$sbatch_dir/$outputPrefix.sbatch.stderr.txt";
	my $sbatch_stdout_path = "$sbatch_dir/$outputPrefix.sbatch.stdout.txt";

	my $src_all_ctss_path = "$baseDir/out/aggregate/$outputPrefix/bed/$outputPrefix.aggregate.collapse.ctss.bed.gz";
	my $src_ung_ctss_path = "$baseDir/out/aggregate/$outputPrefix/bed/$outputPrefix.aggregate.unencoded_G.collapse.ctss.bed.gz";
	my $cp_all_ctss_path = "$all_ctss_dir/$outputPrefix.aggregate.collapse.ctss.bed.gz";
	my $cp_ung_ctss_path = "$ung_ctss_dir/$outputPrefix.aggregate.unencoded_G.collapse.ctss.bed.gz";

	open (SBATCHCMD, ">", $sbatch_sh_path);
	print SBATCHCMD "#!/bin/bash\n";
	print SBATCHCMD "echo \"$runID is started at \$(date)\" >>$start_time_log_path\n";
	print SBATCHCMD "run_script=\"$run_script\"\n";
	print SBATCHCMD "lib_list_path=\"$lib_list_path\"\n";
	print SBATCHCMD "genome=\"$genome\"\n";
	print SBATCHCMD "baseDir=\"$baseDir\"\n";
	print SBATCHCMD "outputPrefix=\"$outputPrefix\"\n";
	print SBATCHCMD 'sh $run_script $lib_list_path $genome $baseDir $outputPrefix'."\n";
	print SBATCHCMD "cp $src_all_ctss_path $cp_all_ctss_path\n";
	print SBATCHCMD "cp $src_ung_ctss_path $cp_ung_ctss_path\n";
	print SBATCHCMD "echo \"$runID is finished at \$(date)\" >>$finish_time_log_path\n";
	close SBATCHCMD;

	print ALLPOOLCTSSLIST join "", (join "\t", ($runID, $cp_all_ctss_path, $cp_ung_ctss_path)), "\n";
	print "$runID sbatch submitted.\n";
	system "sbatch -c 10 -e $sbatch_stderr_path -o $sbatch_stdout_path $sbatch_sh_path";
}
close ALLPOOLCTSSLIST;
