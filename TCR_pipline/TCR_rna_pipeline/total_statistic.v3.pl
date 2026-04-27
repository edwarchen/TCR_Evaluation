#!/usr/local/bin/perl -w
# 去掉flash，直接使用clean_read作为总reads数，merge_rate固定为1，并且修正了clone_reads_ratio的计算方式，确保分母为clean_read。
use strict;
use Cwd;
use Getopt::Long;
my $BEGIN=time();
use Data::Dumper;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);

my $programe_dir=basename($0);
my $path=dirname($0);

# ------------------------------------------------------------------
# GetOptions
# ------------------------------------------------------------------
my($od,$key,$in,$chain);
GetOptions(
    "h|?" => \&help,
    "o:s" => \$od,
    "k:s" => \$key,
    "i:s" => \$in,
    "c:s" => \$chain,
) || &help;

&help unless ($in && $key && $od && $chain);

sub help {
    print <<"Usage End.";
    Description:
    Usage:
        -i          input dir    force 
        -o          output dir   force 
        -k          key of sample force 
        -c          chain (TRA/TRB/BOTH) force
        -h          Help document
Usage End.
    exit;
}

###############Time
my $Time_Start = sub_format_datetime(localtime(time()));
print "\n[$Time_Start] $Script start ... \n\n";

########################################################

# ---------- 处理链类型 ----------
$chain = uc($chain);
if($chain ne "TRA" && $chain ne "TRB" && $chain ne "BOTH"){
    die "Error: -c must be TRA / TRB / BOTH\n";
}

my @chains;
if($chain eq "BOTH"){
    @chains = ("TRA","TRB");
}else{
    @chains = ($chain);
}

# ---------- fastp 文件 ----------
my $fastp_json="$in/Cleanfq/${key}.json";
open (IN,$fastp_json) or die "cannot open $fastp_json \n";

# ---------- 输出 ----------
open (OUT,">$od/${key}.result.index.txt") or die "cannot open output\n";

print OUT "sample\tchain\traw_read\traw_q30\tclean_rate\tclean_read\tclean_q30\tmerge_rate\ttotal_clone_reads\tclone_reads_ratio\tlack_clone_ratio\tClonereads\tSamplingReads\tClone_AA_count\tTop_Freq\tClonality_index\tEvenness\tShannon_Weiner\tSimpson\tMax_Shannon_Weiner\n";

# ---------- fastp 统计 ----------
my $line=0; 
my ($raw_q30,$clean_q30,$raw_read,$clean_read);

while(<IN>){
    chomp;
    $line++;
    if($line==4){ ($raw_read)=(split/:|,/,$_)[1]; }
    if($line==9){ ($raw_q30)=(split/:|,/,$_)[1]; }
    if($line==15){ ($clean_read)=(split/:|,/,$_)[1];}
    if($line==20){ ($clean_q30)=(split/:|,/,$_)[1]; }
}
close(IN);

my $clean_rate=sprintf ("%.4f",$clean_read/$raw_read);

# ===============================
# 🔥 核心修改：不再使用FLASH
# ===============================
my $total_reads = $clean_read;     # MiXCR标准：reads
my $merge_PE_read = $total_reads;  # 兼容旧变量
my $merge_rate = 1;                # 无实际意义，仅占位

########################################################
# 🔥 TRA / TRB 循环
########################################################

foreach my $chain_type (@chains){

    my $convert="$in/Map_Clone_Analysis/convert.${key}.clonotypes.${chain_type}.txt";

    unless(-e $convert){
        warn "Warning: $convert not found, skip $chain_type\n";
        next;
    }

    open (IN1,$convert) or die "cannot open $convert \n";

    my $total_clone_reads=0;
    while(<IN1>){
        chomp;
        next if(/count/);
        my @aa=split/\t/;
        $total_clone_reads+=$aa[0];
    }
    close(IN1);

    # 🔥 核心比例（已修正）
    my $clone_reads_ratio = 0;
    if($total_reads > 0){
        $clone_reads_ratio = sprintf ("%.4f",$total_clone_reads/$total_reads);
    }

    my $lack_clone_ratio = sprintf("%.4f", 1 - $clone_reads_ratio);

    # ---------- 重新读取 ----------
    open (FILE,$convert) or die "cannot open $convert \n";
    <FILE>;

    my %cdr3_aa_new;
    my %info;
    my $sum = 0;
    my $num=1000000;

    while(<FILE>){
        chomp;
        my(@aa)=split/\t/;
        $cdr3_aa_new{$aa[2]}=$aa[0];
        $sum+=$aa[0];
        $info{$aa[2]}=join("\t",@aa[2..10])."\n";
    }
    close(FILE);

    ################## sampling ##################
    my %random; 
    my %random_clone_info; 
    my %random_clone_count;

    if($num >= $sum){
        $random{$_} = 1 for(1..$sum);
    }else{
        while(1){
            my $f = int rand($sum);
            $random{($f+1)}=1;
            last if(scalar keys %random == $num);
        }
    }

    my $samplingreads=$num>=$sum ? $sum : $num;

    my $flag = 1;
    my $max_freq=0;

    for(keys %cdr3_aa_new){
        my $f = 0;
        for(my $i=$flag; $i<$flag+$cdr3_aa_new{$_} ; $i++){
            if(exists $random{$i}){ $f++; }
        }
        $flag += $cdr3_aa_new{$_};

        if($f>0){
            my $freq=$f/$samplingreads;
            if($freq>$max_freq){$max_freq=$freq;}
            $random_clone_count{$_}=$f;
            $random_clone_info{$_}="$f\t$freq\t$info{$_}";
        }
    }

    ################ 多样性 ################
    my $total_count=0;
    my $Simpson_AA_1=0;
    my $Shannon_Weiner_AA=0;
    my %AA;
    my $clonality=0;

    foreach my $key1 (keys %random_clone_count) {
        my @aa=split/\t/,$random_clone_info{$key1};

        next if($aa[0]<1);

        $total_count+=$aa[0];
        $clonality+=$aa[1]*(log($aa[1]));

        if(exists $AA{$aa[3]}){
            $AA{$aa[3]}+=$aa[0];
        }else{
            $AA{$aa[3]}=$aa[0];
        }
    }

    foreach my $key2 (keys %AA) {
        $Simpson_AA_1+=$AA{$key2}*($AA{$key2}-1);
        my $per=$AA{$key2}/$total_count;
        $Shannon_Weiner_AA+=-$per*(log($per));
    }

    my $Count_AA = keys %AA;
    my $max_Shannon_Weiner_AA=log($Count_AA);
    my $Simpson_AA=1-($Simpson_AA_1/($total_count*($total_count-1)));
    my $Evenness= $Count_AA>1 ? $Shannon_Weiner_AA/$max_Shannon_Weiner_AA : 0;
    my $clonality_index=$Count_AA>1 ? 1+$clonality/(log($Count_AA)) : 0;

    # ---------- 输出 ----------
    print OUT "$key\t$chain_type\t$raw_read\t$raw_q30\t$clean_rate\t$clean_read\t$clean_q30\t$merge_rate\t$total_clone_reads\t$clone_reads_ratio\t$lack_clone_ratio\t$sum\t$samplingreads\t$Count_AA\t$max_freq\t$clonality_index\t$Evenness\t$Shannon_Weiner_AA\t$Simpson_AA\t$max_Shannon_Weiner_AA\n";
}

close(OUT);

########################################################
###############Time
my $Time_End = sub_format_datetime(localtime(time()));
print "\n[$Time_End] $Script end, ";
&Runtime($BEGIN);

###########subs

sub Runtime {
    my ($t1)=@_;
    my $t=time()-$t1;
    print "Total elapsed time: ${t}s\n";
}

sub sub_format_datetime {
    my($sec, $min, $hour, $day, $mon, $year) = @_;
    sprintf("%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $day, $hour, $min, $sec);
}