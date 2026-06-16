#!/usr/bin/perl -w
use strict;

my $kallisto_index_hsh_ref = {
    'SALA_table5_partialYes.ENST_chro' => '/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/exosome_sensitivity/kallisto_20240904/index/chr_table5_partial_yes_detected.alone_allNeuron_THP1t5.rm_rRNA',
	'SALA_table5_partialYes.ENST' => '/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/exosome_sensitivity/kallisto_20240904/index/table5.final.partial_yes_detected.alone.rm_rRNA'
};

my $outDir          = "/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/exosome_sensitivity/kallisto_20240904/kallisto_quantification_exosome";
my $in_fq_list_path = "/analysisdata/fantom6/Interactome/ONT.CAGE.satellite/dorado_run/exosome_sensitivity/exo_fastq_path.txt";

system "mkdir -pm 744 $outDir";
system "cp $0 $outDir";

my $start_time_log_path  = "$outDir/00_start.time.log.txt";
my $finish_time_log_path = "$outDir/00_finish.time.log.txt";
system "echo \"========== All runs are started at \$(date)\" ========== >$start_time_log_path\n";
system "echo \"========== All runs are started at \$(date)\" ========== >$finish_time_log_path\n";

my $lib_info_hsh_ref = {};

open my $FQLIST, "<", $in_fq_list_path or die "Cannot open $in_fq_list_path: $!";
while (<$FQLIST>) {
    chomp;
    next if /^#/;
    my ($libID, $fq_path) = split /\t/;
    $lib_info_hsh_ref->{$libID} = $fq_path;
}
close $FQLIST;

my $num_lib = keys %{$lib_info_hsh_ref};
print "$num_lib libraries read\n";

foreach my $indexID (keys %{$kallisto_index_hsh_ref}) {
    my $kallisto_index = $kallisto_index_hsh_ref->{$indexID};
    foreach my $libID (sort keys %{$lib_info_hsh_ref}) {
        my $fq_path = $lib_info_hsh_ref->{$libID};

        die "$fq_path of $libID is not defined\n" if not -s $fq_path;
        my $sub_outDir   = "$outDir/out/$indexID/$libID";
        my $kallistr_dir = "$sub_outDir/";
        system "mkdir -pm 744 $kallistr_dir";

        my $kallisto_cmd = join " ", (
            "kallisto quant",
            "-i $kallisto_index",
            "-o $kallistr_dir",
            "--single",
            "-l 200",
            "-s 50",
            "--rf-stranded",
            "$fq_path"
        );

        my $sbatch_sh_path      = "$sub_outDir/$indexID.$libID.sbatch.cmd.sh";
        my $sbatch_stderr_path  = "$sub_outDir/$indexID.$libID.sbatch.stderr.txt";
        my $sbatch_stdout_path  = "$sub_outDir/$indexID.$libID.sbatch.stdout.txt";

        open my $SBATCHCMD, ">", $sbatch_sh_path or die "Cannot open $sbatch_sh_path for writing: $!";
        print $SBATCHCMD "#!/bin/bash\n";
        print $SBATCHCMD "echo \"$indexID.$libID is started at \$(date)\" >>$start_time_log_path\n";
        print $SBATCHCMD "$kallisto_cmd;\n";
        print $SBATCHCMD "echo \"$indexID.$libID is finished at \$(date)\" >>$finish_time_log_path\n";
        close $SBATCHCMD;

        print "$indexID.$libID sbatch submitted.\n";
        system "sbatch -c 10 -e $sbatch_stderr_path -o $sbatch_stdout_path $sbatch_sh_path";
    }
}