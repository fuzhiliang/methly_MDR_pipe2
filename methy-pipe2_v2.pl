#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Cwd;
use FindBin qw($Bin);
use Getopt::Long;
use lib $Bin;
my $conf_file;
GetOptions
(
	"conf|c:s" => \$conf_file,
);
unless( $conf_file )
{
	usage();
}

our ( $BSAligner, $R );
our ( $bs_index, $list_chr_len, $genome_w_fa, $genome_c_fa, $hg_3mer, $TSS, $info, $hg);
our ( $out_dir, $op, $seq_format, $seq_mode, $used_cycles, $sequenc_tot_cycle, $merge, $bin_size_CpG, $bin_size_nonCpG,
		$mismatch, $min_ins, $max_ins, $hit_mode, $thread, $need_trim, $discard_5_x_cycles,$Phred);

# load and check configuration file
my $conf = load_conf( $conf_file );

warn "Use mode : $seq_mode\n";

my $cwd = cwd();
$out_dir = "$cwd/$out_dir" unless $out_dir=~ /^\//;
$info = "$cwd/$info" unless $info =~ /^\//;

our $perl_prog = "$Bin/perl_prog";
our $cpp_prog  = "$Bin/cpp_prog";
our $util_prog = "$Bin/util";
our $R_prog    = "$Bin/R_prog";

my $pe = 0;
$pe = 1 if $seq_mode eq 'pe';
# load INFO file
print STDERR "Loading info file : $info ...\n";
open INFO, "$info" || die( "Error : open INFO file failed : $!" );
my @all_sample;
my @mdrs_group;
my (%desc, %read1, %read2);
while( <INFO> )
{
	next if /^#/;
	next if /^\s+$/;
	chomp;
	my @l = split /\t/;
	my $mdrs;
	if ($l[0] eq 'MDR'){
		$mdrs="$l[1]\t$l[2]";
		push @mdrs_group,$mdrs;
	}else{
		$l[0] =~ s/[\/\\\s]+/_/g;
		my $key = "$l[0]_$l[1]";
		$desc{ $key }  = $l[2];
		$read1{ $key } = rel2abs(\$l[3]);
		warn "$l[3]\n";
		if( $pe )
		{
			$read2{$key} = rel2abs(\$l[4]);
			warn "$l[4]\n";
		}
		else
		{
			$read2{$key} = '';
		}
		push @all_sample, $key;
	}
}
close INFO;

if( $#all_sample==0 && $merge==1 )
{
	print STDERR "WARNING : You only have one sample BUT you use 'merge' option!\n";
	$merge = 0;
}

# prepare directory and qsub files
prepare( $conf_file );

# prepare makefile

# step 1 : fastq statistics and base composition plot
my $fastq_statistics = '';
my @work;

foreach( @all_sample )
{
	$fastq_statistics .= "$op\_logs/$op.$_.read1.fqstat.finished : $read1{$_}\n";
	$fastq_statistics .= "\t$cpp_prog/fqstatistics $used_cycles $read1{$_} > $op\_logs/$op.$_.read1.fqstat && touch $op\_logs/$op.$_.read1.fqstat.finished\n";
	$fastq_statistics .= "$op\_logs/$op.$_.read1.base.composition.pdf.finished : $op\_logs/$op.$_.read1.fqstat.finished\n";
	$fastq_statistics .= "\t$R --slave --args $op\_logs/$op.$_.read1.fqstat $op\_logs/$op.$_.read1 < $R_prog/plot_base_composition.R && touch $op\_logs/$op.$_.read1.base.composition.pdf.finished\n";
	push @work, "$op\_logs/$op.$_.read1.fqstat.finished $op\_logs/$op.$_.read1.base.composition.pdf.finished";
	if( $pe )
	{
		$fastq_statistics .= "$op\_logs/$op.$_.read2.fqstat.finished : $read2{$_}\n";
		$fastq_statistics .= "\t$cpp_prog/fqstatistics $used_cycles $read2{$_} > $op\_logs/$op.$_.read2.fqstat && touch $op\_logs/$op.$_.read2.fqstat.finished\n";

		$fastq_statistics .= "$op\_logs/$op.$_.read2.base.composition.pdf.finished : $op\_logs/$op.$_.read2.fqstat.finished\n";
		$fastq_statistics .= "\t$R --slave --args $op\_logs/$op.$_.read2.fqstat $op\_logs/$op.$_.read2 < $R_prog/plot_base_composition.R && touch $op\_logs/$op.$_.read2.base.composition.pdf.finished\n";
		push @work, "$op\_logs/$op.$_.read2.fqstat.finished $op\_logs/$op.$_.read2.base.composition.pdf.finished";
	}
}
$fastq_statistics = "$op\_logs/base.composition.finished : ". join(' ', @work) . "\n\ttouch $op\_logs/base.composition.finished\n". $fastq_statistics;

# step 2 : fastq trimming and BSAligner alignment
my $trim_PE="$cpp_prog/trim_adaptor_PE_v2.2";
my $trim_SE="$cpp_prog/trim_adaptor_SE_v3";
my $bsaligner = '';
if( $pe )	# paired-end data
{
	foreach( @all_sample )
	{
		if ($read1{$_}=~/\.gz$/){
			`gunzip $read1{$_} -c > $op\_logs/$op.$_\_1.fq && touch $op\_logs/$op.$_.R1.gunzip.finished\n` unless (-f "$op\_logs/$op.$_.R1.gunzip.finished");
			`gunzip $read2{$_} -c > $op\_logs/$op.$_\_2.fq && touch $op\_logs/$op.$_.R2.gunzip.finished\n` unless (-f "$op\_logs/$op.$_.R1.gunzip.finished");
			$read1{$_}="$op\_logs/$op.$_\_1.fq";
			$read2{$_}="$op\_logs/$op.$_\_2.fq";
		}
		if($need_trim == 1 )
		{
			$bsaligner.= "$op\_logs/$op.$_.trim.finished : $read1{$_} $read2{$_}\n";

			$bsaligner.= "\t$trim_PE $discard_5_x_cycles $used_cycles  $read1{$_} $read2{$_} $op\_logs/$op.$_ 20 $Phred && touch $op\_logs/$op.$_.trim.finished\n";
			
			$bsaligner .= "$op\_logs/$op.$_.bsalign.finished : $op\_logs/$op.$_.trim.finished\n";
			$bsaligner .= "\t$BSAligner -a $op\_logs/$op.$_.R1.fq -b $op\_logs/$op.$_.R2.fq -D $bs_index -r $hit_mode -v $mismatch -m $min_ins -x $max_ins -o $op\_alignment/$op.$_.bsalign -2 $op\_logs/$op.$_.unpaired -p $thread 2>$op\_logs/$op.$_.log && touch $op\_logs/$op.$_.bsalign.finished\n";
		}
		else
		{
			$bsaligner .= "$op\_logs/$op.$_.bsalign.finished : \n";
			$bsaligner .= "\t$BSAligner -a $read1{$_} -b $read2{$_} -D $bs_index -r $hit_mode -v $mismatch -m $min_ins -x $max_ins -o $op\_alignment/$op.$_.bsalign -2 $op\_logs/$op.$_.unpaired -p $thread 2>$op\_logs/$op.$_.log && touch $op\_logs/$op.$_.bsalign.finished\n";
		}
	}
}
else	# single-end data
{
	foreach( @all_sample )
	{
		if($need_trim == 1 )
		{
			$bsaligner.= "$op\_logs/$op.$_.trim.finished : $read1{$_}\n";
			$bsaligner.= "\t$trim_SE $discard_5_x_cycles $used_cycles $read1{$_} $op\_logs/$op.$_ 26 $Phred && touch $op\_logs/$op.$_.trim.finished\n";
			$bsaligner .= "$op\_logs/$op.$_.bsalign.finished : $op\_logs/$op.$_.trim.finished\n";
			$bsaligner .= "\t$BSAligner -a $op\_logs/$op.$_.trimmed.fq -D $bs_index -r $hit_mode -v $mismatch -o $op\_alignment/$op.$_.bsalign -p $thread 2>$op\_logs/$op.$_.log && touch $op\_logs/$op.$_.bsalign.finished\n";
		}
		else
		{
			$bsaligner .= "$op\_logs/$op.$_.bsalign.finished : \n";
			$bsaligner .= "\t$BSAligner -a $read1{$_} -D $bs_index -r $hit_mode -v $mismatch -o $op\_alignment/$op.$_.bsalign -p $thread 2>$op\_logs/$op.$_.log && touch $op\_logs/$op.$_.bsalign.finished\n";			
		}
	}
}

# step 3 : deal with merge option
my $merge_data = '';
my @call_sample;
if( $merge )	# merge data
{
	$merge_data .= "$op\_logs/$op.merged.bsalign.finished :";
	foreach( @all_sample )
	{
		$merge_data .= " $op\_logs/$op.$_.bsalign.finished ";
	}
	if($pe)
	{
		$merge_data .= "\n\t$cpp_prog/removeDuplicate_PE_bsalign $list_chr_len $op\_alignment/$op.merged ";
	}
	else
	{
		$merge_data .= "\n\t$cpp_prog/removeDuplicate_SE_bsalign $list_chr_len $op\_alignment/$op.merged ";
	}
	
	my $bsalign_log="";
	foreach( @all_sample )
	{
		$merge_data .= " $op\_alignment/$op.$_.bsalign ";
		$bsalign_log.=" $op\_logs/$op.$_.log";
	}
	
	$merge_data .= " && perl $perl_prog/stat_bsalign_merged.pl $seq_mode $op\_logs/$bsalign_log >$op\_logs/$op.merged.log && touch $op\_logs/$op.merged.bsalign.finished\n";

	@call_sample = ( 'merged' );
}
else
{
	@call_sample = @all_sample;
}

# step 4 : remove duplicate
my $remove_dup = '';
if( $pe )
{
	foreach( @call_sample )
	{
		unless($merge)
		{
			$remove_dup .= "$op\_logs/$op.$_.rmdup.bsalign.finished : $op\_logs/$op.$_.bsalign.finished\n";
			$remove_dup .= "\t$cpp_prog/removeDuplicate_PE_bsalign $list_chr_len $op\_logs/$op.$_ $op\_alignment/$op.$_.bsalign && mv $op\_logs/$op.$_.W.rmdup.bsalign $op\_logs/$op.$_.C.rmdup.bsalign $op\_alignment && touch $op\_logs/$op.$_.rmdup.bsalign.finished\n";
		}
		else
		{
			$remove_dup .= $merge_data;
		}

		$remove_dup .= "$op\_logs/$op.$_.W.length.hist.pdf.finished : $op\_logs/$op.$_.rmdup.bsalign.finished\n";
		$remove_dup .= "\t$R --slave --args $op\_logs/$op.$_.W.len.hist $op\_logs/$op.$_.W < $R_prog/plot_len_hist.R && touch $op\_logs/$op.$_.W.length.hist.pdf.finished\n";
		$remove_dup .= "$op\_logs/$op.$_.C.length.hist.pdf.finished : $op\_logs/$op.$_.rmdup.bsalign.finished\n";
		$remove_dup .= "\t$R --slave --args $op\_logs/$op.$_.C.len.hist $op\_logs/$op.$_.C < $R_prog/plot_len_hist.R && touch $op\_logs/$op.$_.C.length.hist.pdf.finished\n";
	}
}
else
{
	foreach( @call_sample )
	{
		unless($merge)
		{
			$remove_dup .= "$op\_logs/$op.$_.rmdup.bsalign.finished : $op\_logs/$op.$_.bsalign.finished\n";
			$remove_dup .= "\t$cpp_prog/removeDuplicate_SE_bsalign $list_chr_len $op\_logs/$op.$_ $op\_alignment/$op.$_.bsalign && mv $op\_logs/$op.$_.W.rmdup.bsalign $op\_logs/$op.$_.C.rmdup.bsalign $op\_alignment && touch $op\_logs/$op.$_.rmdup.bsalign.finished\n";
		}
		else
		{
			$remove_dup .= $merge_data;
		}
	}
}

# step 5, 6 : bsalign statistics and methyaltion call
my $meth_call = '';
foreach( @call_sample )
{
	$meth_call .= "$op\_logs/$op.$_.W.call.finished : $op\_logs/$op.$_.rmdup.bsalign.finished\n";
	$meth_call .= "\tmkdir -p $op\_logs/bsalign_statistics/$_ && $cpp_prog/meth_call $seq_mode $genome_w_fa $op\_alignment/$op.$_.W.rmdup.bsalign $op\_logs/bsalign_statistics/$_/$op.$_.W && mv $op\_logs/bsalign_statistics/$_/$op.$_.W.call $op\_meth_call && perl -lane 'print if \$\$F[-1]=~/C:G:/i' $op\_meth_call/$op.$_.W.call > $op\_meth_call/$op.$_.W.CpG.call && touch $op\_logs/$op.$_.W.call.finished\n";

	$meth_call .= "$op\_logs/$op.$_.C.call.finished : $op\_logs/$op.$_.rmdup.bsalign.finished\n";
	$meth_call .= "\tmkdir -p $op\_logs/bsalign_statistics/$_ && $cpp_prog/meth_call $seq_mode $genome_c_fa $op\_alignment/$op.$_.C.rmdup.bsalign $op\_logs/bsalign_statistics/$_/$op.$_.C && mv $op\_logs/bsalign_statistics/$_/$op.$_.C.call $op\_meth_call && touch $op\_logs/$op.$_.C.call.finished\n";

	$meth_call .= "$op\_logs/$op.$_.C.rev.call.finished : $op\_logs/$op.$_.C.call.finished\n";

	$meth_call .= "\t$cpp_prog/reverse_crick_metCall $list_chr_len $op\_meth_call/$op.$_.C.call $op\_meth_call/$op.$_.C.rev.call && perl -lane 'print if \$\$F[-1]=~/C:G:/i'  $op\_meth_call/$op.$_.C.rev.call > $op\_meth_call/$op.$_.C.rev.CpG.call && rm $op\_meth_call/$op.$_.C.call && touch $op\_logs/$op.$_.C.rev.call.finished\n";	
	
	##profile CpG in 1Mb win
	$meth_call .= "$op\_logs/$op.$_.profiling.overall.CpG.density.finished : $op\_logs/$op.$_.W.call.finished $op\_logs/$op.$_.C.rev.call.finished\n";
	$meth_call .="\tperl $perl_prog/profile_mC_mergedStrand.pl 1e6 $op\_meth_call/$op.$_.W.CpG.call $op\_meth_call/$op.$_.C.rev.CpG.call > $op\_meth_density/$op.$_.1Mb.overall.CpG.met.density && touch $op\_logs/$op.$_.profiling.overall.CpG.density.finished\n";
	
	$meth_call .= "$op\_logs/$op.$_.profiling.repeat.CpG.density.finished : $op\_logs/$op.$_.W.call.finished $op\_logs/$op.$_.C.rev.call.finished\n";
    $meth_call .="\tperl $perl_prog/profile_mC_mergedStrand_rep.pl 1e6 $op\_meth_call/$op.$_.W.CpG.call $op\_meth_call/$op.$_.C.rev.CpG.call > $op\_meth_density/$op.$_.1Mb.repeat.CpG.met.density && touch $op\_logs/$op.$_.profiling.repeat.CpG.density.finished\n";
	##met per chr
	$meth_call .= "$op\_logs/$op.$_.met.per.chr.finished : $op\_logs/$op.$_.W.call.finished $op\_logs/$op.$_.C.rev.call.finished\n";
	$meth_call .="\tperl $perl_prog/calc_met_per_chr.pl $op\_meth_call/$op.$_.W.CpG.call $op\_meth_call/$op.$_.C.rev.CpG.call >$op\_logs/$op.$_.met.per.chr && touch $op\_logs/$op.$_.met.per.chr.finished\n";
}

# step 7 : statistics
my $statistics = '';
foreach( @call_sample )
{
	$statistics .= "$op\_logs/$op.$_.W.meth.stat.finished : $op\_logs/$op.$_.W.call.finished\n";
	$statistics .= "\tperl $perl_prog/create_met_stat_report.pl $hg_3mer $op\_meth_call/$op.$_.W.call > $op\_logs/$op.$_.W.meth.stat && touch $op\_logs/$op.$_.W.meth.stat.finished\n";
	$statistics .= "$op\_logs/$op.$_.C.meth.stat.finished : $op\_logs/$op.$_.C.rev.call.finished\n";
	$statistics .= "\tperl $perl_prog/create_met_stat_report.pl $hg_3mer $op\_meth_call/$op.$_.C.rev.call > $op\_logs/$op.$_.C.meth.stat && touch $op\_logs/$op.$_.C.meth.stat.finished\n";

	$statistics .= "$op\_logs/$op.$_.final.report.finished : $op\_logs/$op.$_.rmdup.bsalign.finished $op\_logs/$op.$_.W.meth.stat.finished $op\_logs/$op.$_.C.meth.stat.finished\n";
	$statistics .= "\tperl $perl_prog/bsalign_report_v3.pl $sequenc_tot_cycle $op\_logs/$op.$_.read.length.distr  $list_chr_len $op\_logs/$op.$_.log $op\_logs/$op.$_.W.dup.stat $op\_logs/$op.$_.C.dup.stat >$op\_logs/$op.$_.bsalign.report && $perl_prog/final_report_v2.pl $hg_3mer $op\_logs/$op.$_.W.meth.stat $op\_logs/$op.$_.C.meth.stat $op\_logs/$op.$_.W.dup.stat $op\_logs/$op.$_.C.dup.stat >$op\_logs/$op.$_.final.report 2>$op\_logs/$op.$_.context.distr && touch $op\_logs/$op.$_.final.report.finished\n";
}

# step 8 : plot, methylation distribution
my $plots = '';
foreach( @call_sample )
{
	# plot TSS
	$plots .= "$op\_logs/$op.$_.W.TSS.bin200.CpG.density.finished : $op\_logs/$op.$_.W.call.finished\n";
	$plots .= "\tperl $perl_prog/calc_tss_density.pl CpG $TSS 200 $op\_meth_call/$op.$_.W.call > $op\_logs/$op.$_.W.TSS.bin200.CpG.density && touch $op\_logs/$op.$_.W.TSS.bin200.CpG.density.finished\n";
	$plots .= "$op\_logs/$op.$_.CpG.W.bin200.TSS.pdf.finished : $op\_logs/$op.$_.W.TSS.bin200.CpG.density.finished\n";
	$plots .= "\t$R --slave --args $op\_logs/$op.$_.W.TSS.bin200.CpG.density $op\_logs/$op.$_.CpG.W CpG < $R_prog/plot_TSS.R && touch $op\_logs/$op.$_.CpG.W.bin200.TSS.pdf.finished\n";
		# for C, firstly revserse the positions
	$plots .= "$op\_logs/$op.$_.C.TSS.bin200.CpG.density.finished : $op\_logs/$op.$_.C.rev.call.finished\n";
	$plots .= "\t$perl_prog/calc_tss_density.pl CpG $TSS 200 $op\_meth_call/$op.$_.C.rev.call > $op\_logs/$op.$_.C.TSS.bin200.CpG.density && touch $op\_logs/$op.$_.C.TSS.bin200.CpG.density.finished\n";
	$plots .= "$op\_logs/$op.$_.CpG.C.bin200.TSS.pdf.finished : $op\_logs/$op.$_.C.TSS.bin200.CpG.density.finished\n";
	$plots .= "\t$R --slave --args $op\_logs/$op.$_.C.TSS.bin200.CpG.density $op\_logs/$op.$_.CpG.C CpG < $R_prog/plot_TSS.R && touch $op\_logs/$op.$_.CpG.C.bin200.TSS.pdf.finished\n";

	# calculate density
	$plots .= "$op\_logs/$op.$_.W.density.finished : $op\_logs/$op.$_.W.call.finished\n";
	$plots .= "\tperl $perl_prog/calc_context_density.pl $bin_size_CpG $bin_size_nonCpG $list_chr_len $op\_logs/$op.$_.W $op\_meth_call/$op.$_.W.call && touch $op\_logs/$op.$_.W.density.finished\n";
	$plots .= "$op\_logs/$op.$_.C.density.finished : $op\_logs/$op.$_.C.rev.call.finished\n";
	$plots .= "\t$perl_prog/calc_context_density.pl $bin_size_CpG $bin_size_nonCpG $list_chr_len $op\_logs/$op.$_.C $op\_meth_call/$op.$_.C.rev.call && touch $op\_logs/$op.$_.C.density.finished\n";

	$plots .= "$op\_logs/$op.$_.context.bar.pdf.finished : $op\_logs/$op.$_.final.report.finished\n";
	$plots .= "\t$R --slave --args $op\_logs/$op.$_.context.distr $op\_logs/$op.$_ < $R_prog/plot_bar.R && touch $op\_logs/$op.$_.context.bar.pdf.finished\n";

	$plots .= "$op\_logs/$op.$_.mC.percentage.across.context.pdf.finished : $op\_logs/$op.$_.final.report.finished\n";
	$plots .= "\t$R --slave --args $op\_logs/$op.$_.context.distr $op\_logs/$op.$_ < $R_prog/plot_mC_pct_across_context.R && touch $op\_logs/$op.$_.mC.percentage.across.context.pdf.finished\n";

	$plots .= "$op\_logs/$op.$_.bin$bin_size_CpG.context.CpG.density.pdf.finished : $op\_logs/$op.$_.W.density.finished $op\_logs/$op.$_.C.density.finished\n";
	$plots .= "\t$R --slave --args $op\_logs/$op.$_.W.bin$bin_size_CpG.CpG.density $op\_logs/$op.$_.C.bin$bin_size_CpG.CpG.density CpG $op\_logs/$op.$_.bin$bin_size_CpG < $R_prog/plot_CpG_density_allchr.R && touch $op\_logs/$op.$_.bin$bin_size_CpG.context.CpG.density.pdf.finished\n";
}


# write makefile
open MK, ">$out_dir/makefile" || die( "$!" );

print MK "methylation :";
foreach( @call_sample )
{
	print MK " $op\_logs/$op.$_.finished";
}
print MK "\n\t", html(), "\n";

foreach( @call_sample )
{
	if( $pe )
	{
		print MK "$op\_logs/$op.$_.finished : $op\_logs/$op.$_.C.TSS.bin200.CpG.density.finished $op\_logs/$op.$_.W.TSS.bin200.CpG.density.finished $op\_logs/$op.$_.profiling.overall.CpG.density.finished $op\_logs/$op.$_.profiling.repeat.CpG.density.finished $op\_logs/base.composition.finished $op\_logs/$op.$_.W.length.hist.pdf.finished $op\_logs/$op.$_.C.length.hist.pdf.finished $op\_logs/$op.$_.context.bar.pdf.finished $op\_logs/$op.$_.mC.percentage.across.context.pdf.finished $op\_logs/$op.$_.bin$bin_size_CpG.context.CpG.density.pdf.finished $op\_logs/$op.$_.CpG.W.bin200.TSS.pdf.finished $op\_logs/$op.$_.CpG.C.bin200.TSS.pdf.finished $op\_logs/$op.$_.met.per.chr.finished\n";
	}
	else
	{
		print MK "$op\_logs/$op.$_.finished : $op\_logs/$op.$_.C.TSS.bin200.CpG.density.finished $op\_logs/$op.$_.W.TSS.bin200.CpG.density.finished $op\_logs/$op.$_.profiling.overall.CpG.density.finished $op\_logs/$op.$_.profiling.repeat.CpG.density.finished $op\_logs/base.composition.finished $op\_logs/$op.$_.context.bar.pdf.finished $op\_logs/$op.$_.mC.percentage.across.context.pdf.finished $op\_logs/$op.$_.bin$bin_size_CpG.context.CpG.density.pdf.finished $op\_logs/$op.$_.CpG.W.bin200.TSS.pdf.finished $op\_logs/$op.$_.CpG.C.bin200.TSS.pdf.finished $op\_logs/$op.$_.met.per.chr.finished\n";
	}
	print MK "\ttouch $op\_logs/$op.$_.finished\n";
}
print MK $bsaligner, "\n";
print MK $fastq_statistics, "\n";
#print MK $merge_data, "\n";
print MK $remove_dup, "\n";
print MK $meth_call, "\n";
print MK $statistics, "\n";
print MK $plots, "\n";

close MK;

=cut
print STDERR "\n\nCreate makefile done.\n\n";
print STDERR "********** run mode 1 **********\n",
			 "The output is in $out_dir;\n",
			 "Please run the makefile like:\n\n",
			 "\tmake -j 8 -f makefile\n\n",
			 "********************************\n\n";

print STDERR "********** run mode 2 **********\n",
			 "If you want to run the makefile based on the Sun Grid Engine (SGE) with 2 nodes, please run like:\n\n",
			 "\tsh qsub.sh 2 makefile\n\n",
			 "********************************\n";
=cut
system "cd $out_dir && make -j 8 -f $out_dir/makefile && touch $out_dir/makefile.finished " unless (-f "$out_dir/makefile.finished");

#MDRs
my $mdrs='';
my $ANNO="/home/fuzl/soft/methy_pipe/methy-pipe2/bed_files";
print "$#mdrs_group\t@mdrs_group\n";
if ($mdrs_group[0] && $merge == 0){
	foreach( @call_sample ){
		$mdrs.="cd  $out_dir/$op\_meth_call && sh $Bin/utils/split_meth_call/split_meth_call.sh $op.$_.W.CpG.call $op.$_.C.rev.CpG.call $_\_split $_ \n";
#		`cd  $out_dir/$op\_meth_call && sh $Bin/utils/split_meth_call/split_meth_call.sh $op.$_.W.CpG.call $op.$_.C.rev.CpG.call $_\_split $_ `;
	#sh /home/fuzl/soft/methy_pipe/methy-pipe2/utils/split_meth_call/split_meth_call.sh test.PW396w_7.W.CpG.call test.PW396w_7.C.rev.CpG.call PW396w_split PW396w
	}
	foreach(@mdrs_group){
		my ($control,$treat)=split "\t",$_;
		$mdrs.="cd $out_dir/MDRs && $Bin/utils/DMR_calling/auto_DMR.biomarker.sh $out_dir/$op\_meth_call/$control\_split/$control $out_dir/$op\_meth_call/$treat\_split/$treat $control\_refto_$treat \n";
	}
	$mdrs.="cd $out_dir/MDRs &&  perl $Bin/utils/DMR_anno/dmr_anno.pl $ANNO/iGenome.hg19.revised.gff3 all.hyper.call.filtered >all.hyper.call.filtered.anno2gene.xls \n";
	$mdrs.="cd $out_dir/MDRs &&  perl $Bin/utils/DMR_anno/dmr_anno.pl $ANNO/iGenome.hg19.revised.gff3 all.hypo.call.filtered >all.hypo.call.filtered.anno2gene.xls \n";

	open MDR ,">$out_dir/MDRs.sh";
	print MDR "$mdrs";
	system "sh $out_dir/MDRs.sh";
	close MDR; 
}
print "\nThe project is done\n\n";
####################################subroutine
# prepare directory and qsub.sh, qmake.sh
sub prepare
{
	my $conf_file = shift;
	if( -e $out_dir )
	{
		print STDERR "\n\nWARNING : directory '$out_dir' EXISTS !\n\n";
	}
	else
	{
		system "mkdir -p $out_dir";
	}
	system "cp $conf_file $out_dir/conf";

	unless( -e "$out_dir/$op\_summary" )
	{
		system( "mkdir -p $out_dir/$op\_summary" );
	}
	unless( -e "$out_dir/$op\_summary/image" )
	{
		system( "mkdir -p $out_dir/$op\_summary/image" );
	}

	unless( -e "$out_dir/$op\_logs" )
	{
		system( "mkdir -p $out_dir/$op\_logs" );
	}
	
	unless( -e "$out_dir/$op\_alignment" )
	{
		system( "mkdir -p $out_dir/$op\_alignment" );
	}

	unless( -e "$out_dir/$op\_meth_call" )
	{
		system( "mkdir -p $out_dir/$op\_meth_call" );
	}

	unless( -e "$out_dir/$op\_meth_density" )
	{
		system( "mkdir -p $out_dir/$op\_meth_density" );
	}
	unless( -e "$out_dir/MDRs" )
	{
		system( "mkdir -p $out_dir/MDRs" );
	}	
	open OUT, ">$out_dir/qsub.sh";
	my $qsub='
if [ $# -ne 2 ]
then
	echo "usage: sh qsub.sh <node_num> <makefile>"
	exit 1
fi
node_num=$1
mk=$2
qsub -cwd -v PATH -q all.q -pe make $node_num qmake.sh $mk
';
	print OUT $qsub;
	close OUT;

	open OUT, ">$out_dir/qmake.sh";
	my $qmake='
#!/bin/sh
qmake -inherit -- -f $1
';
	print OUT $qmake;
	close OUT;
}

sub html
{
	my $smy;
	$smy .= "cp $info $out_dir/$op\_summary/info ";
	$smy .= "; cd $out_dir/$op\_summary ";
	$smy .= "; cp $out_dir/$op\_logs/*pdf .; cp $out_dir/$op\_logs//*.png image/ ";
	$smy .= "; cp ../conf conf ";
	$smy .= "; cp $out_dir/$op\_logs/*.report . ";
	$smy .= "; perl $perl_prog/makeHTML_v4.pl $merge > summary.html ";
	$smy .= "; cd ..";
	return $smy;
}

sub compress
{
	my $mk .= "\ncompress:\n";
	$mk .= "\techo Gzip is running at the backgroud to compress the bsalign file, please give it some time. ";
	$mk .= "&& gzip *.bsalign *pileup &\n";

	return $mk;
}

sub clean
{
	my $mk .= "\nclean:\n";
	$mk .= "\trm -f *\n";

	return $mk;
}

sub load_conf
{
	my $conffile = shift;
	my %conf;

	# read configuration file
	open CONF, "$conffile" || die( "$!" );
	while( <CONF> )
	{
		chomp;
		next unless /\S+/;
		next if /^#/;
		warn "$_\n";
		my ($key, $value) = split;
		$conf{uc($key)} = $value;
	}
	close CONF;

	# validate configuration file
	# program
	$BSAligner = rel2abs(\$conf{BSALIGNER}) || conf_uncomplete( 'BSALIGNER' );
	$BSAligner = "$BSAligner/BSAligner" if -d $conf{BSALIGNER};
	$R = $conf{R} || conf_uncomplete( 'R' );
	$R = "$conf{R}/R" if( -d $conf{R} );

	# file location
	$bs_index = rel2abs(\$conf{BS_INDEX}) || conf_uncomplete( 'BS_INDEX' );
	$list_chr_len = rel2abs(\$conf{LIST_CHR_LEN}) || conf_uncomplete( 'LIST_CHR_LEN' );
	$genome_w_fa = rel2abs(\$conf{GENOME_W_FA}) || conf_uncomplete( 'GENOME_W_FA' );
	$genome_c_fa = rel2abs(\$conf{GENOME_C_FA}) || conf_uncomplete( 'GENOME_C_FA' );
	$hg_3mer = rel2abs(\$conf{HG_3MER}) || conf_uncomplete( 'HG_3MER' );
	$TSS = rel2abs(\$conf{TSS}) || conf_uncomplete( 'TSS' );
	$info = rel2abs(\$conf{INFO}) || conf_uncomplete( 'INFO' );
	#$hg = exists $conf{HG} ? $conf{HG} || 19;
	$discard_5_x_cycles= $conf{DISCARD_5_X_CYCLES} || 0;
	# options
	$out_dir = $conf{OUT_DIR} || '.';
	$op = $conf{OUT_PREFIX} || 'Methylation';
	$seq_format = lc( $conf{SEQ_FORMAT} || 'fq' );
	$seq_format = 'fq' if $seq_format ne 'fa';
	$seq_mode = lc( $conf{SEQ_MODE} || 'pe' );
	$seq_mode = 'pe' unless $seq_mode eq 'se';
	$used_cycles = $conf{USED_CYCLES} || 50;
	$sequenc_tot_cycle = $conf{SEQUENC_TOT_CYCLE} || 100;
	$merge = $conf{MERGE} || $conf{MERGE_BAM} || 0;
	$bin_size_CpG = $conf{BIN_SIZE_CPG} || 100e3;
	$bin_size_nonCpG = $conf{BIN_SIZE_NONCPG} || 100e3;
	$mismatch = $conf{MISMATCH} || 2;
	$min_ins = $conf{MIN_INS} || 0;
	$max_ins = $conf{MAX_INS} || 600;
	$hit_mode = $conf{HIT_MODE} || 0;
	$thread = $conf{THREAD} || 1;
	$need_trim = exists $conf{TRIM} ? $conf{TRIM} : 1;
	$Phred= $conf{Phred} || 38;
	
	return \%conf;
}

sub conf_uncomplete
{
	my $msg = shift;
	warn( "Configuration file uncompelte : missing key $msg.\n" );
	exit 1;
}

sub rel2abs
{
	my $path_ref=shift;
	return '' unless defined $$path_ref;

	my $current_dir=cwd();
	$$path_ref = "$current_dir/$$path_ref" unless $$path_ref=~/^\//;
	return $$path_ref;
}

sub usage
{
	my $u=qq(
                        ===== METHYLATION ANALYSIS PIPELINE v2.01 =====
contact: jiangpeiyong\@cuhk.edu.hk
2012-3-28 23:00:00

Note1* this pipeline use Burrows-Wheeler Transform (BWT) to design BS-Seq aligner, named BSAligner.
Usage: perl $0 --conf conf

-c --conf                    specify the conf file
);
	warn "$u\n";
	exit 1;
}


